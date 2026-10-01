# data_visualization
# 
# Data exploration
# 

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(readr)
# library(cmdstanr)
# library(rstanarm)
# library(posterior)
# library(bayesplot)
# library(sf)
# library(rnaturalearth)

source("SWFSC_Stranding_fcns.R")

# The database was accessed using SQL_data_extraction.R on the following date.
# All tables were stored as .csv files - this decision might not have been the
# best because there are some fields that had commas in their entries. 
data.extraction.date <- "2026-10-01"

# 
# table.names <- c("_Animal", "_Morphology", "_Age", "_Reproduction",
#                  "_Bone", "Code_Maturity")

tables. <- read.tables(data.extraction.date)

tables.$Morph.2 %>%
  filter(Genus.f == "Tursiops") %>%
  mutate(Species.f = as.factor(Species)) %>% 
  left_join(tables.$Animal %>%
              select(Specimen, Sex, Month, Day), by = "Specimen") -> table.Morph.2.Tursiops

write.csv(table.Morph.2.Tursiops,
          file = "data//Tursiops_stranding.csv")
