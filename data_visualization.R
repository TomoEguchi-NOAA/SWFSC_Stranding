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

MCMC.params <- list(n.samples = 25000,
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

# All column types are stored in the SWFSC_Stranding_fcns.R script

# There were some extensive editing (mostly deleting entries) for the Species
# table. So, I'm not going to change the data file from 2026-05-12
table.Species <- read_csv(file = paste0("Data/tblSpecies.csv"),
                          col_types = species.col.types) %>%
  rename(SpeciesID = SpCode) %>%
  select(-c(EditDate, EditUser, RecordCreationDate))

table.Morph <- read_csv(file = paste0("Data/tbl_Morphology_", 
                                      data.extraction.date, ".csv"),
                        col_types = morph.col.types) %>%
  select(-c(EditDate, EditUser, RecordCreationDate))

table.Animal <- read_csv(file = paste0("Data/tbl_Animal_", 
                                       data.extraction.date, ".csv"),
                         col_types = animal.col.types) %>%
  select(-c(EditDate, EditUser, RecordCreationDate))


table.Age <- read_csv(file = paste0("Data/tbl_Age_", 
                                    data.extraction.date, ".csv"),
                      col_types = age.col.types) %>%
  select(-c(EditDate, EditUser, RecordCreationDate))

table.Age %>%
  filter(IsAnalysisQuality == "Y") %>%
  left_join(table.Animal, by = "Specimen") -> table.Age.1

table.Repro <- read_csv(file = paste0("Data/tbl_Reproduction_", 
                                      data.extraction.date, ".csv"),
                        col_types = repro.col.types) %>%
  select(-c(EditDate, EditUser, RecordCreationDate))

table.Weight <- read_csv(file = paste0("Data/tbl_Weight_", 
                                       data.extraction.date, ".csv"),
                         col_types = weight.col.types) %>%
  select(-c(EditDate, EditUser, RecordCreationDate))

# Remove a fin whale that had very small Lab measurement
# The filter removes if ratio = NA. Need to explicitly keep NAs. 
table.Morph  %>%
  select(Specimen, TotalLength_LAB, TotalLength_FIELD, IsStandardTL_FIELD,
         IsStandardTL_LAB, GIRTHMAX) %>%
  mutate(ratio = TotalLength_FIELD/TotalLength_LAB) %>%
  filter(ratio < 2 | is.na(ratio)) -> table.Morph.1

Field.Lab.lm <- lm(TotalLength_LAB ~ TotalLength_FIELD, data = table.Morph.1)

table.Morph.1 %>%
  left_join(table.Animal %>% 
              select(Specimen, SpeciesID, Year, Latitude, Longitude),
            by = "Specimen") %>% 
  left_join(table.Species %>% 
              select(SpeciesID, Genus, Species, CommonName, SpName), 
            by = "SpeciesID") %>% 
  filter(Genus != "Unid") %>%
  mutate(Genus.f = as.factor(Genus)) %>% 
  left_join(table.Weight %>%
              select(Specimen, Carcass_Intact, Heart, Kidney_R, Kidney_L, Liver), 
            by = "Specimen") %>%
  left_join(table.Age.1 %>%
              select(Specimen, Age, IsAnalysisQuality, EstimationMethod), 
            by = "Specimen") %>%
  left_join(table.Repro %>% 
              select(Specimen, IsMature, MaturityID, IsLactating, IsPregnant,
                     CA_LEFT, CA_RIGHT, TotalCorpora), 
            by = "Specimen") -> table.Morph.2

table.Morph.2 %>%
  filter(Genus.f == "Stenella") %>%
  mutate(Species.f = as.factor(Species))-> table.Morph.2.Stenella

# ggplot(table.Morph.2.Stenella) +
#   geom_point(aes(x = Age, y = TotalLength_FIELD, 
#                  color = Species.f, size = Latitude),
#              alpha = 0.5)

jags.out.filename.age <- "RData/jags_out_Stenella_age_logistic.rds"

if (!file.exists(jags.out.filename.age)){
  
  # simple age-sexual maturity logistic regression:
  table.Morph.2.Stenella %>%
    select(Age, IsMature, Species.f, TotalLength_FIELD) %>%
    na.omit() %>%
    droplevels() %>%
    mutate(IsMature.int = as.integer(IsMature),
           Species.int = as.integer(Species.f)) %>%
    mutate(Mature.idx = ifelse(IsMature.int == 1, 1, 0)) -> Stenella.maturity
  
  sp.names <- paste("S.", levels(Stenella.maturity$Species.f))
  
  jags.params <- c("B", "Sigma", "rho", "mu_beta", "loglik")
  # Define an identity matrix for the Wishart scale matrix
  R_matrix <- matrix(c(1, 0, 0, 1), nrow = 2)
  
  jags.data <- list(
    X = Stenella.maturity$Age,          
    mature = Stenella.maturity$Mature.idx,    
    species_idx = Stenella.maturity$Species.int, 
    N = nrow(Stenella.maturity),
    n.sp = max(Stenella.maturity$Species.int),
    R = R_matrix
  )
  
  jags.model <- "models/model_logistic_regression.jags"
  
  tic <- Sys.time()
  jm.Stenella.age <- jagsUI::jags(jags.data,
                                  inits = NULL,
                                  parameters.to.save= jags.params,
                                  model.file = jags.model,
                                  n.chains = MCMC.params$n.chains,
                                  n.burnin = MCMC.params$n.burnin,
                                  n.thin = MCMC.params$n.thin,
                                  n.iter = MCMC.params$n.samples,
                                  DIC = T,
                                  parallel=T)
  
  toc <- Sys.time()
  Rhat.Stenella.age <- rank.normalized.R.hat(jm.Stenella.age$samples,
                                             params = "^mu_beta|^B|^Sigma|^rho",
                                             MCMC.params = MCMC.params)
  
  jm.logistic.summary.age <- summary.logistic(jags.out = jm.Stenella.age,
                                              jags.data = jags.data, 
                                              xvar = "Age",
                                              sp.names = sp.names)

  jags.out.age <- list(jags.model = jags.model,
                       jags.data = jags.data,
                       MCMC.params = MCMC.params,
                       jags.out = jm.Stenella.age,
                       Rhat = Rhat.Stenella.age,
                       logistic.summary = jm.logistic.summary.age,
                       Run.Date = Sys.Date(),
                       Run.Time = toc - tic,
                       System = Sys.getenv())
  
  saveRDS(jags.out.age, 
          file = jags.out.filename.age)
} else {
  jags.out.age <- read_rds(jags.out.filename.age)
} 

jags.out.filename.length <- "RData/jags_out_Stenella_length_logistic.rds"

if (!file.exists(jags.out.filename.length)){
  
  # Do a similar analysis with length:
  jags.data <- list(X = Stenella.maturity$TotalLength_FIELD,          
                    mature = Stenella.maturity$Mature.idx,    
                    species_idx = Stenella.maturity$Species.int, 
                    N = nrow(Stenella.maturity),
                    n.sp = max(Stenella.maturity$Species.int),
                    R = R_matrix)
  
  jags.model <- "models/model_logistic_regression.jags"
  
  jm.Stenella.length <- jagsUI::jags(jags.data,
                                     inits = NULL,
                                     parameters.to.save= jags.params,
                                     model.file = jags.model,
                                     n.chains = MCMC.params$n.chains,
                                     n.burnin = MCMC.params$n.burnin,
                                     n.thin = MCMC.params$n.thin,
                                     n.iter = MCMC.params$n.samples,
                                     DIC = T,
                                     parallel=T)
  
  Rhat.Stenella.length <- rank.normalized.R.hat(jm.Stenella.length$samples,
                                                params = "^mu_beta|^B|^Sigma|^rho",
                                                MCMC.params = MCMC.params)
  
  jm.logistic.summary.length <- summary.logistic(jags.out = jm.Stenella.length,
                                                 jags.data = jags.data, 
                                                 xvar = "Length",
                                                 sp.names = sp.names)
  
  jags.out.length <- list(jags.model = jags.model,
                          jags.data = jags.data,
                          MCMC.params = MCMC.params,
                          jags.out = jm.Stenella.length,
                          Rhat = Rhat.Stenella.length,
                          logistic.summary = jm.logistic.summary.length,
                          Run.Date = Sys.Date(),
                          Run.Time = toc - tic,
                          System = Sys.getenv())
  
  saveRDS(jags.out.length, 
          file = jags.out.filename.length)
} else {
  jags.out.length <- read_rds(jags.out.filename.length)
} 

table.Morph.2 %>%
    filter(Genus.f == "Delphinus") %>%
    mutate(Species.f = as.factor(Species))-> table.Morph.2.Delphinus
  
  # ggplot(table.Morph.2.Delphinus) +
  #   geom_point(aes(x = Age, y = TotalLength_FIELD, 
  #                  color = Species.f, size = Latitude),
  #              alpha = 0.5)
  
  table.Morph.2 %>%
    filter(Genus.f == "Tursiops") %>%
  mutate(Species.f = as.factor(Species))-> table.Morph.2.Tursiops

# ggplot(table.Morph.2.Tursiops) +
#   geom_point(aes(x = Age, y = TotalLength_FIELD, 
#                  color = Latitude))
