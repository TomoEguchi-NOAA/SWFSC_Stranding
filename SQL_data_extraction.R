# Access stranding and life history databases and extract morphological and
# life history data. This project is for a Hollings scholar (Raija Hammond)
# in Summer 2026.
# 
rm(list = ls())
library(RODBC)
library(tidyverse)
library(readr)
library(ggplot2)

source("SWFSC_Stranding_fcns.R")

# get Species table from Common 
Common.con <- connection.string("SWFSCCommon") 
Common  <- odbcDriverConnect(Common.con)
#on.exit(odbcClose(Common))  # close connection if R crashes

Common.tables <- sqlTables(Common) %>%
  filter(TABLE_TYPE == "TABLE")

Common %>%
  sqlQuery('select * from tblSpecies') %>% #-> tmp
  #filter(SubOrder == "CETACEA") %>%
  select(ID, SubOrder, Genus, Species, CommonName,
         SpName, SpCode) %>%
  mutate(SpeciesID = SpCode) -> tbl.Species

odbcClose(Common)

write.csv(tbl.Species,
          file = paste0("Data/table_species_", Sys.Date(), ".csv"))

# Get the lifehistory database
MMLH.con <- connection.string("MMLH") 
MMLH.2019 <- odbcDriverConnect(MMLH.con)
#on.exit(odbcClose(MMLH.2019))  # close connection if R crashes
MMLH.info <- odbcGetInfo(MMLH.2019)
MMLH.tables <- sqlTables(MMLH.2019)

MMLH.tables %>%
  filter(TABLE_TYPE == "TABLE") -> MMLH.table.names

table.names <- c("Animal", "Morphology", "Age", "Reproduction",
                 "Teeth", "Gonads", "Bone", "Code_Maturity",
                 "InvOsteology")

table.list <- list()

for (k in 1:length(table.names)){
  table.name <- paste0("tbl_", table.names[k])
  MMLH.2019 %>%
    sqlQuery(paste0("select * from ", table.name)) -> table.list[[k]]
  
  # write.csv(table.list[[k]],
  #           file = paste0("Data/", table.name, "_", Sys.Date(), ".csv"))
  
}

names(table.list) <- table.names

# The following will not work - need to pull out by names 2026-05-11
table.list[[grep("Animal", table.names)]] %>%
  select(Specimen, Year, Month, Day, Latitude, Latitude_Precision,
         Latitude_Precision_Unit, Longitude, Longitude_Precision,
         Longitude_Precision_Unit, SpeciesID, Sex) %>%
  left_join(tbl.Species, by = "SpeciesID")-> tbl.Animal

# Select Delphinus and Tursiops
tbl.Animal %>%
  filter(Genus == "Delphinus" |
         Genus == "Tursiops") -> tbl.Animal.dolphins

table.list[[grep("Morphology", table.names)]] %>%
  select(Specimen, IsStandardTL_LAB, TotalLength_LAB, 
         IsStandardTL_FIELD, TotalLength_FIELD, STOANUS,
         STOGENSLIT, STOUMBIL, STOTHRGROO, STODOFINTIP,
         STOANTDOR, STOFLIPPER, STOEAR, STOEYE, STOGAPE,
         STOBLOHOLE, STOMELAPEX, GIRTHMAX) -> tbl.Morphology

tbl.Animal.dolphins %>%
  left_join(tbl.Morphology, by = "Specimen") -> tbl.dolphins.morphology

table.list[[grep("Age", table.names)]] %>%
  select(-c(Comments, EditDate, EditUser, RecordCreationDate)) -> tbl.Age

tbl.dolphins.morphology %>%
  left_join(tbl.Age, by = "Specimen") -> tbl.dolphins.morphology.age

tbl.dolphins.morphology.age %>%
  filter(IsStandardTL_LAB == "Y" | IsStandardTL_FIELD == "Y") -> tbl.dolphins.StandardTL

tbl.dolphins.morphology.age %>%
  filter(IsStandardTL_LAB == "N" | IsStandardTL_FIELD == "N") -> tbl.dolphins.NoStandardTL

# Create a linear model to predict the STL from other measurements:
# Doesn't work because so many missing data... 
# lm.STL.1 <- lm(TotalLength_FIELD ~ STOANUS + STOGENSLIT + STOUMBIL + STOTHRGROO + STODOFINTIP + STOANTDOR + STOFLIPPER + STOEAR + STOEYE + STOGAPE + STOBLOHOLE + STOMELAPEX + GIRTHMAX,
#                data = tbl.dolphins.StandardTL)

odbcClose(MMLH.2019)

Tissue.con <- connection.string("TissueArchive")
Tissue  <- odbcDriverConnect(Tissue.con)

#on.exit(odbcClose(Tissue))  # close connection if R crashes
odbcClose(Tissue)


# in case there are stray databases open:
odbcCloseAll()


