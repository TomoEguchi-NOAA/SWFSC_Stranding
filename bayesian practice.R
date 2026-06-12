
rm(list=ls())

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
  #filter(Genus.f == "Stenella") %>%
  mutate(Species.f = as.factor(Species)) %>% 
  left_join(tables.$Animal %>%
              select(Specimen, Sex), by = "Specimen") %>%
  select(Specimen, Genus.f, Species.f, TotalLength_FIELD, TotalLength_LAB) %>%
  rename(FIELD = TotalLength_FIELD,
         LAB = TotalLength_LAB) %>%
  na.omit()-> table.Morph.2

write.csv(table.Morph.2,
          file = "data/TotalLength.csv")

ggplot(table.Morph.2) +
  geom_point(aes(x = LAB, y = FIELD))

jags.data <- list(x = table.Morph.2$LAB,
                  y = table.Morph.2$FIELD,
                  N.samples = nrow(table.Morph.2))

jm. <- jagsUI::jags(jags.data,
                    inits = NULL,
                    parameters.to.save = c("beta0", "beta1", "s"),
                    model.file = "models/model_linear_regression.jags",
                    n.chains = 3,
                    n.burnin = 500,
                    n.thin = 5,
                    n.iter = 5000,
                    DIC = T,
                    parallel=T)


lm.1 <- lm(FIELD ~ LAB, data = table.Morph.2)

tables.$Morph.2 %>%
  filter(Genus.f == "Stenella") %>%
  mutate(Species.f = as.factor(Species)) %>% 
  left_join(tables.$Animal %>%
              select(Specimen, Sex), by = "Specimen") -> table.Morph.2.Stenella.1

# remove a couple of "outliers," which are probably data entry errors:
table.Morph.2.Stenella.1 %>% 
  #filter(Longitude < -78) %>%
  select(TotalLength_FIELD, TotalLength_LAB) %>%
  na.omit()-> table.length



ggplot(table.length) +
  geom_point(aes(x = TotalLength_LAB, y = TotalLength_FIELD))

jags.data <- list(x = table.length$TotalLength_LAB,
                  y = table.length$TotalLength_FIELD,
                  N.samples = nrow(table.length))


jm. <- jagsUI::jags(jags.data,
                    inits = NULL,
                    parameters.to.save = c("beta0", "beta1"),
                    model.file = "models/model_linear_regression.jags",
                    n.chains = 3,
                    n.burnin = 500,
                    n.thin = 5,
                    n.iter = 5000,
                    DIC = T,
                    parallel=T)

lm.1 <- lm(y ~ x, data = jags.data)


x <- seq(-50, 50)
y <- dnorm(x, mean = 0, sd = 10)
plot(x,y)
