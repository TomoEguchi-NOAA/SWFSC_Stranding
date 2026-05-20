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
               select(Specimen, Sex), by = "Specimen") -> table.Morph.2.Stenella

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
                         n.sp = max(Stenella.maturity$Species.int),
                         R = R_matrix)

jags.out.age <- jags.logistic(jags.data = jags.data.age,
                              MCMC.params = MCMC.params.1,
                              out.filename = "RData/jags_out_Stenella_age_logistic.rds",
                              jags.model = "models/model_logistic_regression.jags",
                              jags.params = c("B", "Sigma", "rho", "mu_beta", "loglik"),
                              xvar = "Age",
                              yvar = "Probability of maturity",
                              sp.names)


## The average length at sexual maturity
## 
jags.data.length <- list(X = Stenella.maturity$TotalLength_FIELD,          
                         y = Stenella.maturity$Mature.idx,    
                         species_idx = Stenella.maturity$Species.int, 
                         N = nrow(Stenella.maturity),
                         n.sp = max(Stenella.maturity$Species.int),
                         R = R_matrix)

jags.out.length <- jags.logistic(jags.data = jags.data.length,
                                 MCMC.params = MCMC.params.1,
                                 out.filename = "RData/jags_out_Stenella_length_logistic.rds",
                                 jags.model = "models/model_logistic_regression.jags",
                                 jags.params = c("B", "Sigma", "rho", "mu_beta", "loglik"),
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

jags.data.LAB <- list(X = Length.at.Birth.data$Length,          
                      y = Length.at.Birth.data$Status,    
                      species_idx = Length.at.Birth.data$Species.int, 
                      N = nrow(Length.at.Birth.data),
                      n.sp = max(Length.at.Birth.data$Species.int),
                      R = R_matrix)

jags.out.LAB <- jags.logistic(jags.data = jags.data.LAB,
                              MCMC.params = MCMC.params.2,
                              out.filename = "RData/jags_out_Stenella_LAB_logistic.rds",
                              jags.model = "models/model_logistic_regression.jags",
                              jags.params = c("B", "Sigma", "rho", "mu_beta", "loglik"),
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
         Sex_idx = as.numeric(Sex)) -> Length.at.Age.data

# Define the identity matrix for the Wishart prior
#R_matrix <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)

generate_inits_1sex <- function() {
  list(
    # Logistic parameters
    mu_L0 = runif(1, min = 80, max = 85),
    L0 = runif(3, min = 80, max = 85),
    mu_B1 = runif(1, min = 0.05, max = 0.15),
    B1 = runif(3, min = 0.05, max = 0.15),
    
    # Growth parameters
    mu_tc = runif(1, min = 5, max = 7),
    tc = runif(3, min = 5, max = 7),
    
    mu_a1 = rnorm(1, mean = 0, sd = 0.1), tau_a1 = runif(1, 0.5, 1.5), 
    a1 = runif(3, min = 0.1, max = 0.5),
    
    mu_alpha1 = rnorm(1, mean = 0, sd = 0.1), tau_alpha1 = runif(1, 0.5, 1.5), 
    alpha1 = runif(3, min = 0.1, max = 0.5),
    
    mu_a2 = rnorm(1, mean = 0, sd = 0.1), tau_a2 = runif(1, 0.5, 1.5), 
    a2 = runif(3, min = 0.1, max = 0.5),
    
    mu_alpha2 = rnorm(1, mean = 0, sd = 0.1), tau_alpha2 = runif(1, 0.5, 1.5), 
    alpha2 = runif(3, min = 0.1, max = 0.5),
    
    # Variance parameter
    tau = runif(1, min = 0.5, max = 2.0)
  )
}

jags.out.Laird <- jags.Laird.growth(LAB.data = Length.at.Birth.data,
                                    length.age.data = Length.at.Age.data,
                                    MCMC.params = MCMC.params.2,
                                    out.filename = "RData/jags_out_Stenella_Laird.rds",
                                    jags.model = "models/model_two_phase_Laird_Growth.jags",
                                    jags.params = c("L0", "B1", "a1", "alpha1", "a2", 
                                                    "alpha2", "tc", "sigma", "loglik"),
                                    inits.fcn = generate_inits_1sex)

## Add sex as another factor with additional parameters:

# 1. Create a function to generate sensible starting values
generate_inits_sex <- function() {
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
    
    # Variance
    tau = runif(1, min = 0.5, max = 2.0)
  )
}

jags.out.Laird.sex <- jags.Laird.growth(LAB.data = Length.at.Birth.data,
                                        length.age.data = Length.at.Age.data,
                                        MCMC.params = MCMC.params.2,
                                        out.filename = "RData/jags_out_Stenella_Laird_sex.rds",
                                        jags.model = "models/model_two_phase_Laird_Growth_sex.jags",
                                        jags.params = c("B1", "L0", "s_a1", "s_a2", "s_alpha1", "s_alpha2",
                                                        "s_tc", "a1", "alpha1", "a2", "alpha2", "tc", "loglik"),
                                        inits.fcn = generate_inits_sex)


# Plot estimates:

# 1. Extract the MCMC samples into a usable matrix
# Assuming your coda.samples output is named 'mcmc_samples'
post_samples <- as.matrix(jags.out.age.length$jags.out$samples)

# 2. Create a smooth sequence of ages for the x-axis
age_seq <- seq(0, max(growth_data$age, na.rm = TRUE), length.out = 100)

# 3. Create an empty dataframe to store our plotting data
plot_data <- data.frame()

# 4. Loop through 3 species and 2 sexes
for (s in 1:3) {
  for (x in 1:2) {
    
    # Extract shared Phase 1 and changepoint parameters for species 's'
    L0_samp     <- post_samples[, paste0("L0[", s, "]")]
    a1_samp     <- post_samples[, paste0("a1[", s, "]")]
    alpha1_samp <- post_samples[, paste0("alpha1[", s, "]")]
    tc_samp     <- post_samples[, paste0("tc[", s, "]")]
    
    # Extract sex-specific Phase 2 parameters for species 's' and sex 'x'
    # The extraction string looks like "a2[1,1]", "a2[1,2]", etc.
    a2_samp     <- post_samples[, paste0("a2[", s, ",", x, "]")]
    alpha2_samp <- post_samples[, paste0("alpha2[", s, ",", x, "]")]
    
    # Matrix to hold predictions (Rows = MCMC iterations, Columns = age points)
    pred_matrix <- matrix(NA, nrow = nrow(post_samples), ncol = length(age_seq))
    
    for (i in 1:length(age_seq)) {
      current_age <- age_seq[i]
      
      # Phase 1 calculation
      mu1 <- L0_samp * exp((a1_samp / alpha1_samp) * (1 - exp(-alpha1_samp * current_age)))
      
      # Phase 2 calculation
      Lc <- L0_samp * exp((a1_samp / alpha1_samp) * (1 - exp(-alpha1_samp * tc_samp)))
      mu2 <- Lc * exp((a2_samp / alpha2_samp) * (1 - exp(-alpha2_samp * (current_age - tc_samp))))
      
      # Determine phase
      pred_matrix[, i] <- ifelse(current_age < tc_samp, mu1, mu2)
    }
    
    # 5. Calculate median and 95% credible intervals
    species_sex_df <- data.frame(
      species = s,
      sex = x,
      age = age_seq,
      fit = apply(pred_matrix, 2, quantile, probs = 0.500),
      lwr = apply(pred_matrix, 2, quantile, probs = 0.025),
      upr = apply(pred_matrix, 2, quantile, probs = 0.975)
    )
    
    plot_data <- rbind(plot_data, species_sex_df)
  }
}

# 6. Format the labels for ggplot2
# Predicted curve labels
plot_data$species_label <- factor(plot_data$species, 
                                  levels = 1:3, 
                                  labels = c("S. attenuata", "S. coeruleoalba", "S. longirostris"))
plot_data$sex_label <- factor(plot_data$sex, 
                              levels = c(1, 2), 
                              labels = c("Female", "Male"))

# Raw data labels (assuming 'growth_data$sex' is already 1 for Female, 2 for Male)
growth_data$species_label <- factor(as.numeric(as.factor(growth_data$species)), 
                                    levels = 1:3, 
                                    labels = c("S. attenuata", "S. coeruleoalba", "S. longirostris"))
growth_data$sex_label <- factor(as.numeric(as.factor(growth_data$sex)), 
                                levels = c(1, 2), 
                                labels = c("Female", "Male"))

# 7. Plot using ggplot2
ggplot() +
  # 95% credible interval ribbons mapped to sex
  geom_ribbon(data = plot_data, aes(x = age, ymin = lwr, ymax = upr, fill = sex_label), 
              alpha = 0.3) +
  
  # Median fitted curves mapped to sex
  geom_line(data = plot_data, aes(x = age, y = fit, color = sex_label), 
            linewidth = 1) +
  
  # Raw data points mapped to sex
  geom_point(data = growth_data, aes(x = age, y = length, color = sex_label), 
             alpha = 0.4, size = 1.5) +
  
  # Facet by species
  facet_wrap(~species_label, scales = "free_x") +
  
  # Custom colors to easily distinguish males and females
  scale_color_manual(values = c("Female" = "darkorange", "Male" = "steelblue")) +
  scale_fill_manual(values = c("Female" = "darkorange", "Male" = "steelblue")) +
  
  # Formatting
  theme_minimal() +
  labs(
    title = "Sex-Specific Two-Phase Laird Growth Curves",
    x = "Age (Years / Layers)",
    y = "Length (cm)",
    color = "Sex",
    fill = "Sex",
    subtitle = "Solid lines indicate median fits; shaded regions indicate 95% credible intervals"
  ) +
  theme(
    strip.text = element_text(face = "italic", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )