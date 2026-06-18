# data_visualization
# 
# Data exploration
# 

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(readr)
library(cmdstanr)
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

# 1. Prepare the data
stan_data_mature_age <- list(
  N = nrow(Stenella.maturity),
  N_species = length(sp.names),
  X = Stenella.maturity$Age,
  y = Stenella.maturity$Mature.idx,
  species = Stenella.maturity$Species.int
)

# 2. Compile and Fit the Model
logistic_model <- cmdstan_model("models/model_logistic_regression.stan")

if (!file.exists("RData/mature_age_logistic.rds")){
  mature.age_fit <- logistic_model$sample(
    data = stan_data_mature_age,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000
  )
  
  mature.age_fit$save_object(file = "RData/mature_age_logistic.rds")
  
} else {
  mature.age_fit <- readRDS(file = "RData/mature_age_logistic.rds")
}

if (!file.exists("RData/mature_length_logistic.rds")){
  stan_data_mature_length <- list(
    N = nrow(Stenella.maturity),
    N_species = length(sp.names),
    X = Stenella.maturity$TotalLength_FIELD,
    y = Stenella.maturity$Mature.idx,
    species = Stenella.maturity$Species.int
  )
  
  mature.length_fit <- logistic_model$sample(
    data = stan_data_mature_length,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000
  )
  
  mature.length_fit$save_object(file = "RData/mature_length_logistic.rds")
  
} else {
  mature.length_fit <- readRDS(file = "RData/mature_length_logistic.rds")
}

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

stan_data_length_at_birth <- list(
  N = nrow(Length.at.Birth.data),
  N_species = length(sp.names),
  X = Length.at.Birth.data$Length,
  y = Length.at.Birth.data$Status,
  species = Length.at.Birth.data$Species.int
)

if (!file.exists("RData/length_at_birth_logistic.rds")){
  length.at.birth_fit <- logistic_model$sample(
    data = stan_data_length_at_birth,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000
  )
  
  length.at.birth_fit$save_object(file = "RData/length_at_birth_logistic.rds")
  
} else {
  length.at.birth_fit <- readRDS(file = "RData/length_at_birth_logistic.rds")
}

# 3. View the L0 Estimates
L0_summary <- length.at.birth_fit$summary(variables = "L0")
print(L0_summary)

# 4. Extract the Vectors for Your Growth Model!
# You can now inject these directly into the data list for your main model
mu_L0_data <- L0_summary$mean
sd_L0_data <- L0_summary$sd

table.Morph.2.Stenella %>%
  select(Species.f, Age, TotalLength_FIELD, Sex) %>%
  filter(Species.f != "frontalis") %>%
  filter(Sex != "U") %>%
  droplevels() %>%
  na.omit() %>%
  mutate(species_idx = as.numeric(Species.f),
         Sex_idx = as.numeric(Sex),
         Length = TotalLength_FIELD) -> Length.at.Age.data

# There are some Age = 0 individuals with length > 100. I filter them out because
# they can't be right.
Length.at.Age.data %>% 
  filter(Age == 0 & Length > 100) -> Too.Big.at.Birth

Length.at.Age.data %>% 
  filter(!(Age == 0 & Length > 100)) -> Length.at.Age.data

stan_data_growth <- list(
  N_sexes = 2,
  N_growth = nrow(Length.at.Age.data),
  age_growth = Length.at.Age.data$Age,
  length_growth = Length.at.Age.data$Length,
  sex_growth = Length.at.Age.data$Sex_idx,
  sp_growth = Length.at.Age.data$species_idx,
  N_species = 3,
  has_phase2 = c(1, 0, 1),
  
  # Injecting our posterior summaries as fixed data
  mu_L0_data = mu_L0_data, 
  sd_L0_data = sd_L0_data
)

# Compile and Fit the Model
length_at_age_model <- cmdstan_model("models/model_growth_model_delta.stan")

if (!file.exists("RData/length_at_age_Lairds.rds")){
  length.at.age_fit <- length_at_age_model$sample(
    data = stan_data_growth,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000
  )
  
  length.at.age_fit$save_object(file = "RData/length_at_age_Lairds.rds")
  
} else {
  length.at.age_fit <- readRDS(file = "RData/length_at_age_Lairds.rds")
}
