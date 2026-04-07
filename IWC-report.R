# IWC annual report on Life History database
# 
# 

rm(list = ls())
library(RODBC)
library(tidyverse)
library(readr)

source("SWFSC_Stranding_fcns.R")

save.files <- F

# get Species table from Common 
# The function "connection.string" is in SWFSC_Stranding_fcns.R. Modify it
# when SQL databases move to another location.
Common.con <- connection.string("SWFSCCommon") 
Common  <- odbcDriverConnect(Common.con)
#on.exit(odbcClose(Common))  # close connection if R crashes

# The following line returns an error when "Source"d the script
# but it runs fine when executed by one line.
Common.tables <- sqlTables(Common) %>%
  filter(TABLE_TYPE == "TABLE")

Common %>%
  sqlQuery('select * from tblSpecies') %>% #-> tmp
  #filter(SubOrder == "CETACEA") %>%
  select(ID, SubOrder, Genus, Species, CommonName,
         SpName, SpCode) %>%
  mutate(SpeciesID = SpCode) -> tbl.Species

odbcClose(Common)

# Get the lifehistory database
MMLH.con <- connection.string("MMLH") 
MMLH.2019 <- odbcDriverConnect(MMLH.con)
#on.exit(odbcClose(MMLH.2019))  # close connection if R crashes
MMLH.info <- odbcGetInfo(MMLH.2019)
MMLH.tables <- sqlTables(MMLH.2019)

MMLH.tables %>%
  filter(TABLE_TYPE == "TABLE") -> MMLH.table.names

MMLH.2019 %>%
  sqlQuery('select * from tbl_Animal') %>%  #-> MMLH.2019.all
  select(Specimen, Year, Month, Day, Latitude, Latitude_Precision,
         Latitude_Precision_Unit, Longitude, Longitude_Precision,
         Longitude_Precision_Unit, SpeciesID) %>%
  left_join(tbl.Species, by = "SpeciesID")-> tbl.Animal

odbcClose(MMLH.2019)

# Need # species, # records, Start and end years
n.records <- nrow(tbl.Animal)

# Count the number of entries per suborder
tbl.Animal %>%
  group_by(SubOrder) %>%
  summarise(n = n()) -> n.Order

# Interesting that Sirenia is not entered... it
# returns NA for manatee
tbl.Animal %>%
  filter(is.na(SubOrder)) %>%
  filter(!is.na(CommonName)) -> No.Order

tbl.Animal %>%
  filter(SubOrder == "CETACEA") %>%
  filter(SpName != "Unid") %>%
  filter(!str_detect(SpName, "spp")) %>%
  filter(!str_detect(tolower(CommonName), "unidentified")) %>%
  group_by(SpeciesID) %>%
  summarise(n = n(),
            SpName = first(SpName),
            CommonName = first(CommonName)) %>%
  arrange(desc(n)) -> n.Cetacea.species 

tbl.Animal %>%
  filter(SubOrder == "PINNIPEDIA") %>%
  filter(SpName != "Unid") %>%
  filter(!str_detect(SpName, "spp")) %>%
  filter(!str_detect(tolower(CommonName), "unidentified")) %>%
  group_by(SpeciesID) %>%
  summarise(n = n(),
            SpName = first(SpName),
            CommonName = first(CommonName)) %>%
  arrange(desc(n)) -> n.Pinnipedia.species 


start.year <- min(tbl.Animal$Year, na.rm = T)
end.year <- max(tbl.Animal$Year, na.rm = T)

if (save.files){
  write.csv(n.Cetacea.species, 
            file = paste0("IWC_Cetacea_", date(now()), ".csv"))
  write.csv(n.Pinnipedia.species, 
            file = paste0("IWC_Pinnipedia_", date(now()), ".csv"))
  
}

# close all odbc connections in case there are strays 
odbcCloseAll()
