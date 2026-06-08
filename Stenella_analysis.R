# data_visualization
# 
# Data exploration
# 

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(readr)
library(jagsUI)
library(rstanarm)
library(posterior)
library(bayesplot)
library(sf)
library(rnaturalearth)

source("SWFSC_Stranding_fcns.R")
options(mc.cores = parallel::detectCores())

MCMC.params.1 <- list(n.samples = 30000,
                      n.thin = 100,
                      n.burnin = 10000,
                      n.chains = 5)

# The database was accessed using SQL_data_extraction.R on the following date.
# All tables were stored as .csv files - this decision might not have been the
# best because there are some fields that had commas in their entries. 
data.extraction.date <- "2026-05-11"

# 
# table.names <- c("_Animal", "_Morphology", "_Age", "_Reproduction",
#                  "_Bone", "Code_Maturity")

# Define an identity matrix for the Wishart scale matrix
R_matrix <- matrix(c(1, 0, 0, 1), nrow = 2)

tables. <- read.tables(data.extraction.date)

tables.$Morph.2 %>%
  filter(Genus.f == "Stenella") %>%
  mutate(Species.f = as.factor(Species)) %>% 
  left_join(tables.$Animal %>%
              select(Specimen, Sex), by = "Specimen") -> table.Morph.2.Stenella.1

# remove a couple of "outliers," which are probably data entry errors:
table.Morph.2.Stenella.1 %>% 
  filter(Longitude < -78) -> table.Morph.2.Stenella.2


# 1. Convert your dataframe into a spatial 'sf' object
# Replace 'porpoise_data' with your actual dataframe name, 
# and ensure 'lon' and 'lat' match your column names.
# CRS 4326 is the standard coordinate reference system for GPS lat/lon.
points_sf <- st_as_sf(table.Morph.2.Stenella.2, 
                      coords = c("Longitude", "Latitude"), crs = 4326)

# 2. Download the high-resolution global coastline
# We use scale = "large" for the most accurate coastal boundaries
coastline <- ne_coastline(scale = "large", returnclass = "sf")

# 3. Check which points are within 25 km (25,000 meters) of the coast
# Because the CRS is 4326, the 'sf' package uses spherical geometry (s2) 
# and automatically calculates the distance in meters.
within_25km <- st_is_within_distance(points_sf, coastline, dist = 25000)

# 4. Create your binary covariate
# 'within_25km' returns a list. If a point touches the coast, its list length is > 0.
table.Morph.2.Stenella.2$is_coastal <- ifelse(lengths(within_25km) > 0, 1, 0)

# View how many dolphins were flagged as coastal!
table(table.Morph.2.Stenella.2$is_coastal)

# Remove coastal (only 58 data points)
table.Morph.2.Stenella.2 %>%
  filter(is_coastal == 0) -> table.Morph.2.Stenella

# simple age-sexual maturity logistic regression:
table.Morph.2.Stenella %>%
  select(Age, IsMature, Species.f, TotalLength_FIELD, IsStandardTL_FIELD, 
         Year, Latitude, Longitude) %>%
  na.omit() %>%
  droplevels() %>%
  mutate(IsMature.int = as.integer(IsMature),
         Species.int = as.integer(Species.f)) %>%
  mutate(Mature.idx = ifelse(IsMature.int == 1, 1, 0)) -> Stenella.maturity

sp.names <- paste("S.", levels(Stenella.maturity$Species.f))

## The average age at sexual maturity (A50)
jags.data.age <- list(X = Stenella.maturity$Age,          
                      y = Stenella.maturity$Mature.idx,    
                      species_idx = Stenella.maturity$Species.int, 
                      N = nrow(Stenella.maturity),
                      n.sp = max(Stenella.maturity$Species.int))

# 1. Update the function to accept a starting guess for the 50% mark
generate_inits_logistic <- function(start_p50) {
  list(
    # Center the random guesses around whatever start_p50 you provide
    mu_L0 = runif(1, min = start_p50 - 10, max = start_p50 + 10),
    L0 = runif(3, min = start_p50 - 10, max = start_p50 + 10),
    
    # Keep the slope small and positive
    mu_slope = runif(1, min = 0.01, max = 0.1),
    slope = runif(3, min = 0.01, max = 0.1)
  )
}

n_chains <- MCMC.params.1$n.chains
inits.age <- lapply(1:n_chains, function(i) generate_inits_logistic(start_p50 = 10))

jags.out.age <- jags.logistic(jags.data = jags.data.age,
                              MCMC.params = MCMC.params.1,
                              out.filename = "RData/jags_out_Stenella_age_logistic.rds",
                              jags.model = "models/model_logistic_regression.jags",
                              jags.params = c("L0", "slope", "mu_L0", "mu_slope", "tau_L0", "tau_slope"),
                              inits = inits.age,
                              xvar = "Age",
                              yvar = "Probability of maturity",
                              sp.names)


## The average length at sexual maturity
## 
jags.data.length <- list(X = Stenella.maturity$TotalLength_FIELD,          
                         y = Stenella.maturity$Mature.idx,    
                         species_idx = Stenella.maturity$Species.int, 
                         N = nrow(Stenella.maturity),
                         n.sp = max(Stenella.maturity$Species.int))

n_chains <- MCMC.params.1$n.chains
inits.length <- lapply(1:n_chains, function(i) generate_inits_logistic(start_p50 = 188))

jags.out.length <- jags.logistic(jags.data = jags.data.length,
                                 MCMC.params = MCMC.params.1,
                                 out.filename = "RData/jags_out_Stenella_length_logistic.rds",
                                 jags.model = "models/model_logistic_regression.jags",
                                 jags.params = c("L0", "slope", "mu_L0", "mu_slope", 
                                                 "tau_L0", "tau_slope"),
                                 inits = inits.length,
                                 xvar = "Length",
                                 yvar = "Probability of maturity",
                                 sp.names)

## Age and Growth
# Age at birth
# 1. Minimum length by species
table.Morph.2.Stenella %>%
  select(Species.f, TotalLength_FIELD) %>%
  group_by(Species.f) %>%
  na.omit() %>%  
  summarize(length = min(TotalLength_FIELD, na.rm = T),
            n = n()) -> min.length

# 2. Maximum length of fetus by species
table.Morph.2.Stenella %>%
  select(Species.f, FetusLength_Standard) %>%
  group_by(Species.f) %>%
  na.omit() %>%
  summarize(length = max(FetusLength_Standard, na.rm = T),
            n = n()) -> max.length.fetus

# Combine the non-fetus and fetus data and run a logistic regression:
table.Morph.2.Stenella %>%
  select(Species.f, TotalLength_FIELD) %>%
  filter(Species.f != "frontalis") %>%
  droplevels() %>%
  na.omit() %>%
  mutate(Species.int = as.integer(Species.f),
         Status = 1) %>%
  rename(Length = TotalLength_FIELD) %>%
  filter(Length < 125) -> postnatal

table.Morph.2.Stenella %>%
  select(Species.f, FetusLength_Standard) %>%
  filter(Species.f != "frontalis") %>%
  droplevels() %>%
  na.omit() %>%
  mutate(Species.int = as.integer(Species.f),
         Status = 0) %>%
  rename(Length = FetusLength_Standard) -> prenatal

# No fetus record for frontalis
Length.at.Birth.data <- rbind(postnatal, prenatal) 

# Convergence issues. Increased the number of samples
MCMC.params.2 <- list(n.samples = 120000,
                      n.thin = 100,
                      n.burnin = 80000,
                      n.chains = 5)

# Test MCMC setup for debugging
# MCMC.params.2 <- list(n.samples = 1200,
#                       n.thin = 2,
#                       n.burnin = 800,
#                       n.chains = 5)

jags.data.LAB <- list(X = Length.at.Birth.data$Length,          
                      y = Length.at.Birth.data$Status,    
                      species_idx = Length.at.Birth.data$Species.int, 
                      N = nrow(Length.at.Birth.data),
                      n.sp = max(Length.at.Birth.data$Species.int))

n_chains <- MCMC.params.2$n.chains
inits.birth <- lapply(1:n_chains, function(i) generate_inits_logistic(start_p50 = 82.5))

jags.out.LAB <- jags.logistic(jags.data = jags.data.LAB,
                              MCMC.params = MCMC.params.2,
                              out.filename = "RData/jags_out_Stenella_LAB_logistic.rds",
                              jags.model = "models/model_logistic_regression.jags",
                              jags.params = c("L0", "slope", "mu_L0", "mu_slope", 
                                              "tau_L0", "tau_slope"),
                              inits = inits.birth,
                              xvar = "Length at Birth",
                              yvar = "Proportion Postnatal",
                              sp.names)

table.Morph.2.Stenella %>%
  select(Species.f, Age, TotalLength_FIELD, Sex) %>%
  filter(Species.f != "frontalis") %>%
  filter(Sex != "U") %>%
  droplevels() %>%
  na.omit() %>%
  mutate(species_idx = as.numeric(Species.f),
         Sex_idx = as.numeric(Sex),
         Length = TotalLength_FIELD) -> Length.at.Age.data

# Define the identity matrix for the Wishart prior
#R_matrix <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

generate_inits_1sex <- function() {
  list(
    # We can still start L0 near our empirical means to help it initialize smoothly
    L0 = mu_L0_est, 
    
    # Growth parameters
    mu_tc = runif(1, min = 5, max = 7),
    tc = runif(3, min = 5, max = 7),
    
    mu_a1 = rnorm(1, mean = 0, sd = 0.1), 
    tau_a1 = runif(1, 0.5, 1.5), 
    a1 = runif(3, min = 0.1, max = 0.5),
    
    mu_alpha1 = rnorm(1, mean = 0, sd = 0.1), 
    tau_alpha1 = runif(1, 0.5, 1.5), 
    alpha1 = runif(3, min = 0.1, max = 0.5),
    
    mu_a2 = rnorm(1, mean = 0, sd = 0.1), 
    au_a2 = runif(1, 0.5, 1.5), 
    a2 = runif(3, min = 0.1, max = 0.5),
    
    mu_alpha2 = rnorm(1, mean = 0, sd = 0.1), 
    tau_alpha2 = runif(1, 0.5, 1.5), 
    alpha2 = runif(3, min = 0.1, max = 0.5),
    
    tau = runif(1, min = 0.5, max = 2.0)
  )
}

# Just growth, without estimating L0 at the same time
post_matrix <- as.matrix(jags.out.LAB$jags.out$samples)

# 1. Calculate the mean L0 for each species
mu_L0_est <- c(
  mean(post_matrix[, "L0[1]"]),
  mean(post_matrix[, "L0[2]"]),
  mean(post_matrix[, "L0[3]"])
)

# 2. Calculate the standard deviation, then convert to precision (tau = 1 / variance)
sd_L0_est <- c(
  sd(post_matrix[, "L0[1]"]),
  sd(post_matrix[, "L0[2]"]),
  sd(post_matrix[, "L0[3]"])
)
tau_L0_est <- 1 / (sd_L0_est^2)

jags.data.growth <- list(
  N_growth = nrow(Length.at.Age.data),
  age_growth = Length.at.Age.data$Age,
  length_growth = Length.at.Age.data$Length,
  sp_growth = Length.at.Age.data$species_idx,
  N_species = 3,
  
  # Injecting our posterior summaries as fixed data
  mu_L0_data = mu_L0_est, 
  tau_L0_data = tau_L0_est
)

# Create the list of lists for 5 chains
# inits_list_growth <- lapply(1:MCMC.params.2$n.chains, 
#                                function(i) generate_inits_growth())

jags.out.Laird <- jags.Laird.growth(jags.data = jags.data.growth,
                                    MCMC.params = MCMC.params.2,
                                    out.filename = "RData/jags_out_Stenella_Laird.rds",
                                    jags.model = "models/model_two_phase_Laird.jags",
                                    jags.params = c("L0", "a1", "alpha1", "a2", "alpha2", "tc",
                                                    "s_a1", "s_a2", "s_alpha1", "s_alpha2",
                                                    "s_tc", "sigma", "loglik"),
                                    inits.fcn = generate_inits_1sex)

## Add sex as another factor with additional parameters:

generate_inits_sex <- function() {
  list(
    L0 = mu_L0_est, 
    
    tc = matrix(runif(6, min = 5, max = 7),
                nrow = 3, ncol = 2),
    
    mu_tc = runif(2, min = 5, max = 7),
    
    a1 = matrix(runif(3, min = 0.1, max = 0.5),
                nrow = 3, ncol = 2),
    alpha1 = matrix(runif(3, min = 0.1, max = 0.5),
                    nrow = 3, ncol = 2),
    
    mu_a2 = matrix(runif(6, min = -0.5, max = 0.5), 
                nrow = 3, ncol = 2),
    
    # Initialize L_inf around the 200cm mark
    L_inf = matrix(runif(6, min = 190, max = 210), 
                   nrow = 3, ncol = 2),
    
    tau = runif(1, min = 0.5, max = 2.0)
  )
}


jags.data.growth.sex <- jags.data.growth
jags.data.growth.sex$sex_growth <- Length.at.Age.data$Sex_idx
jags.data.growth.sex$N_sex <- max(Length.at.Age.data$Sex_idx)
jags.data.growth.sex$has_phase2 = c(1, 0, 1)  # Sp1 = Yes, Sp2 = No, Sp3 = Yes

jags.out.Laird.sex <- jags.Laird.growth(jags.data = jags.data.growth.sex,
                                        MCMC.params = MCMC.params.2,
                                        out.filename = "RData/jags_out_Stenella_Laird_sex_Linf.rds",
                                        jags.model = "models/model_two_phase_Laird_sex_Linf.jags",
                                        jags.params = c("L0", "Lc", "a1", "alpha1", "a2", "alpha2", "tc",
                                                        "s_a1", "s_a2", "s_alpha1", "s_alpha2",
                                                        "s_tc", "sigma", "L_inf", "s_L_inf",
                                                        "mu_a1", "mu_tc", "mu_alpha1", "mu_a2", 
                                                        "loglik"),
                                        inits.fcn = generate_inits_sex)

#Test MCMC setup for debugging
# MCMC.params.2 <- list(n.samples = 1200,
#                       n.thin = 2,
#                       n.burnin = 800,
#                       n.chains = 5)


generate_inits_sex_scaled <- function() {
  list(
    L0 = mu_L0_est, 
    
    tc = matrix(runif(6, min = 0.5, max = 0.7),
                nrow = 3, ncol = 2),
    
    mu_tc = runif(2, min = 0.5, max = 0.7),
    
    a1 = matrix(runif(3, min = 0.1, max = 0.5),
                nrow = 3, ncol = 2),
    alpha1 = matrix(runif(3, min = 0.1, max = 0.5),
                    nrow = 3, ncol = 2),
    
    mu_a2 = matrix(runif(6, min = -0.5, max = 0.5), 
                   nrow = 3, ncol = 2),
    
    # Initialize L_inf around the 200cm mark
    L_inf = matrix(runif(6, min = 190, max = 210), 
                   nrow = 3, ncol = 2),
    
    tau = runif(1, min = 0.5, max = 2.0)
  )
}


## Scaled age growth model
jags.data.growth.sex$age_scaled <- jags.data.growth.sex$age_growth/10

jags.out.Laird.sex.scaled <- jags.Laird.growth(jags.data = jags.data.growth.sex,
                                               MCMC.params = MCMC.params.2,
                                               out.filename = "RData/jags_out_Stenella_Laird_sex_Linf_scaled.rds",
                                               jags.model = "models/model_two_phase_Laird_sex_Linf_scaled.jags",
                                               jags.params = c("L0", "Lc", "a1", "alpha1", "a2", "alpha2", "tc",
                                                               "s_a1", "s_a2", "s_alpha1", "s_alpha2",
                                                               "s_tc", "sigma", "L_inf", "s_L_inf",
                                                               "mu_a1", "mu_tc", "mu_alpha1", "mu_a2", 
                                                               "loglik"),
                                               inits.fcn = generate_inits_sex_scaled)





generate_inits_1sex_2vars <- function() {
  list(
    # We can still start L0 near our empirical means to help it initialize smoothly
    L0 = mu_L0_est, 
    
    # Growth parameters
    mu_tc = runif(1, min = 5, max = 7),
    tc = runif(3, min = 5, max = 7),
    
    mu_a1 = rnorm(1, mean = 0, sd = 0.1), 
    tau_a1 = runif(1, 0.5, 1.5), 
    a1 = runif(3, min = 0.1, max = 0.5),
    
    mu_alpha1 = rnorm(1, mean = 0, sd = 0.1), 
    tau_alpha1 = runif(1, 0.5, 1.5), 
    alpha1 = runif(3, min = 0.1, max = 0.5),
    
    mu_a2 = rnorm(1, mean = 0, sd = 0.1), 
    tau_a2 = runif(1, 0.5, 1.5), 
    a2 = runif(3, min = 0.1, max = 0.5),
    
    mu_alpha2 = rnorm(1, mean = 0, sd = 0.1), 
    tau_alpha2 = runif(1, 0.5, 1.5), 
    alpha2 = runif(3, min = 0.1, max = 0.5),
    
    # Differentiated starting precisions
    tau_juv = runif(1, min = 1.0, max = 3.0),   # Tighter spread
    tau_adult = runif(1, min = 0.1, max = 0.5)  # Wider spread
  )
}

jags.out.Laird.1sex.2vars <- jags.Laird.growth(jags.data = jags.data.growth.sex,
                                               MCMC.params = MCMC.params.2,
                                               out.filename = "RData/jags_out_Stenella_Laird_2vars.rds",
                                               jags.model = "models/model_two_phase_Laird_2vars.jags",
                                               jags.params = c("L0", "a1", "alpha1", "a2", "alpha2", "tc",
                                                               "s_a1", "s_a2", "s_alpha1", "s_alpha2",
                                                               "s_tc", "s_juv", "s_adult", "loglik"),
                                               inits.fcn = generate_inits_1sex_2vars)


generate_inits_sex_2vars <- function() {
  list(
    # Logistic parameters
    mu_L0 = runif(1, min = 80, max = 85),
    L0 = runif(3, min = 80, max = 85),
    mu_B1 = runif(1, min = 0.05, max = 0.15),
    B1 = runif(3, min = 0.05, max = 0.15),
    
    # Phase 1 parameters (Vectors of length 3)
    mu_tc = runif(1, min = 5, max = 7),
    tc = runif(3, min = 5, max = 7),
    mu_a1 = rnorm(1, mean = 0, sd = 0.1), 
    tau_a1 = runif(1, 0.5, 1.5), 
    a1 = runif(3, min = 0.1, max = 0.5),
    mu_alpha1 = rnorm(1, mean = 0, sd = 0.1), 
    tau_alpha1 = runif(1, 0.5, 1.5), 
    alpha1 = runif(3, min = 0.1, max = 0.5),
    
    # Phase 2 parameters (Matrices: 3 rows for species, 2 columns for sex)
    mu_a2 = rnorm(2, mean = 0, sd = 0.1), 
    tau_a2 = runif(1, 0.5, 1.5), 
    a2 = matrix(runif(6, min = 0.1, max = 0.5), nrow = 3, ncol = 2),
    
    mu_alpha2 = rnorm(2, mean = 0, sd = 0.1), 
    tau_alpha2 = runif(1, 0.5, 1.5), 
    alpha2 = matrix(runif(6, min = 0.1, max = 0.5), nrow = 3, ncol = 2),
    
    # Differentiated starting precisions
    tau_juv = runif(1, min = 1.0, max = 3.0),   # Tighter spread
    tau_adult = runif(1, min = 0.1, max = 0.5)  # Wider spread
  )
}


jags.out.Laird.2sex.2vars <- jags.Laird.growth(jags.data = jags.data.growth.sex,
                                               MCMC.params = MCMC.params.2,
                                               out.filename = "RData/jags_out_Stenella_Laird_sex_2vars.rds",
                                               jags.model = "models/model_two_phase_Laird_sex_2vars.jags",
                                               jags.params = c("L0", "a1", "alpha1", "a2", "alpha2", "tc",
                                                               "s_a1", "s_a2", "s_alpha1", "s_alpha2",
                                                               "s_tc", "s_juv", "s_adult", "loglik"),
                                               inits.fcn = generate_inits_sex_2vars)



## Mixture models with offshore-coastal separation have been removed

#Test MCMC setup for debugging
# MCMC.params.2 <- list(n.samples = 1200,
#                       n.thin = 2,
#                       n.burnin = 800,
#                       n.chains = 5)
# 
# generate_inits_mixture <- function() {
#   list(
#     L0 = mu_L0_est, 
#     
#     mu_tc = runif(1, min = 5, max = 7),
#     tc = runif(3, min = 5, max = 7),
#     
#     mu_a1 = rnorm(1, mean = 0, sd = 0.1), 
#     tau_a1 = runif(1, 0.5, 1.5), 
#     a1 = runif(3, min = 0.1, max = 0.5),
#     
#     mu_alpha1 = rnorm(1, mean = 0, sd = 0.1), 
#     tau_alpha1 = runif(1, 0.5, 1.5), 
#     alpha1 = runif(3, min = 0.1, max = 0.5),
#     
#     mu_a2 = rnorm(1, mean = 0, sd = 0.1), 
#     tau_a2 = runif(1, 0.5, 1.5), 
#     a2 = matrix(runif(6, min = 0.1, max = 0.5), nrow = 3, ncol = 2),
#     
#     # Initialize L_inf inside our bounded constraints
#     L_inf = matrix(c(200, 200, 200, 230, 200, 200), nrow = 3, ncol = 2),
#     
#     p_coastal = runif(1, 0.1, 0.9),
#     
#     s_juv = runif(1, min = 1.0, max = 3.0),   
#     s_adult = runif(1, min = 0.1, max = 0.5)  
#   )
# }
# 
# jags.out.Laird.2vars.mix <- jags.Laird.growth(jags.data = jags.data.growth.sex,
#                                                MCMC.params = MCMC.params.2,
#                                                out.filename = "RData/jags_out_Stenella_Laird_2vars_mix.rds",
#                                                jags.model = "models/model_two_phase_Laird_2vars_mix.jags",
#                                                jags.params = c("L0", "a1", "alpha1", "a2", "tc",
#                                                                "s_a1", "s_a2", "s_alpha1", "s_tc",  
#                                                                "s_juv", "s_adult", "z_aux", "p_coastal",
#                                                                "loglik"),
#                                                inits.fcn = generate_inits_mixture)
# 
# 
# generate_inits_mix_sex <- function() {
#   # Build a 3D array for L_inf starting values
#   L_inf_start <- array(200, dim = c(3, 2, 2))
#   
#   # Give S. attenuata sensible starts based on our bounds
#   L_inf_start[1, 1, 1] <- 190  # Offshore Female
#   L_inf_start[1, 1, 2] <- 200  # Offshore Male
#   L_inf_start[1, 2, 1] <- 220  # Coastal Female
#   L_inf_start[1, 2, 2] <- 230  # Coastal Male
#   
#   list(
#     L0 = mu_L0_est, 
#     tc = runif(3, min = 5, max = 7),
#     a1 = runif(3, min = 0.1, max = 0.5),
#     alpha1 = runif(3, min = 0.1, max = 0.5),
#     a2 = matrix(runif(6, min = 0.1, max = 0.5), nrow = 3, ncol = 2),
#     L_inf = L_inf_start,
#     p_coastal = runif(1, 0.1, 0.9),
#     s_juv = runif(1, min = 1.0, max = 3.0),   
#     s_adult = runif(1, min = 0.1, max = 0.5)  
#   )
# }
# 
# jags.out.Laird.sex.2vars.mix <- jags.Laird.growth(jags.data = jags.data.growth.sex,
#                                                   MCMC.params = MCMC.params.2,
#                                                   out.filename = "RData/jags_out_Stenella_Laird_sex_2vars_mix.rds",
#                                                   jags.model = "models/model_two_phase_Laird_sex_2vars_mix.jags",
#                                                   jags.params = c("L0", "a1", "alpha1", "a2", "tc",
#                                                                   "s_a1", "s_a2", "s_alpha1", "s_tc",  
#                                                                   "s_juv", "s_adult", "z_aux", "p_coastal",
#                                                                   "loglik"),
#                                                   inits.fcn = generate_inits_mix_sex)
# 


# 
# # 6. Format the labels for ggplot2
# # Predicted curve labels
# plot_data$species_label <- factor(plot_data$species, 
#                                   levels = 1:3, 
#                                   labels = c("S. attenuata", "S. coeruleoalba", "S. longirostris"))
# plot_data$sex_label <- factor(plot_data$sex, 
#                               levels = c(1, 2), 
#                               labels = c("Female", "Male"))
# 
# # Raw data labels (assuming 'growth_data$sex' is already 1 for Female, 2 for Male)
# growth_data$species_label <- factor(as.numeric(as.factor(growth_data$species)), 
#                                     levels = 1:3, 
#                                     labels = c("S. attenuata", "S. coeruleoalba", "S. longirostris"))
# growth_data$sex_label <- factor(as.numeric(as.factor(growth_data$sex)), 
#                                 levels = c(1, 2), 
#                                 labels = c("Female", "Male"))
# 
# # 7. Plot using ggplot2
# ggplot() +
#   # 95% credible interval ribbons mapped to sex
#   geom_ribbon(data = plot_data, aes(x = age, ymin = lwr, ymax = upr, fill = sex_label), 
#               alpha = 0.3) +
#   
#   # Median fitted curves mapped to sex
#   geom_line(data = plot_data, aes(x = age, y = fit, color = sex_label), 
#             linewidth = 1) +
#   
#   # Raw data points mapped to sex
#   geom_point(data = growth_data, aes(x = age, y = length, color = sex_label), 
#              alpha = 0.4, size = 1.5) +
#   
#   # Facet by species
#   facet_wrap(~species_label, scales = "free_x") +
#   
#   # Custom colors to easily distinguish males and females
#   scale_color_manual(values = c("Female" = "darkorange", "Male" = "steelblue")) +
#   scale_fill_manual(values = c("Female" = "darkorange", "Male" = "steelblue")) +
#   
#   # Formatting
#   theme_minimal() +
#   labs(
#     title = "Sex-Specific Two-Phase Laird Growth Curves",
#     x = "Age (Years / Layers)",
#     y = "Length (cm)",
#     color = "Sex",
#     fill = "Sex",
#     subtitle = "Solid lines indicate median fits; shaded regions indicate 95% credible intervals"
#   ) +
#   theme(
#     strip.text = element_text(face = "italic", size = 12),
#     plot.title = element_text(face = "bold", size = 14),
#     legend.position = "bottom"
#   )