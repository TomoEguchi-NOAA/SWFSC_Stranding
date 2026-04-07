# DA_2025_summary
# Creates summary of the DA event in 2025.

rm(list = ls())
library(tidyverse)
library(readxl)

dat.cols <- cols("Specimen" = col_character(),
                 "DATE" = col_date(format = "%m/%d/%Y"),
                 "Species" = col_character(),
                 "T24hr_Report_Filed" = col_character(),
                 "Level_A" = col_character(),
                 "Level_A_to_Alaina" = col_character(),
                 "Entered_into_tblPhoneLog" = col_character(),
                 "HI_entered_into_national_db" = col_character(),
                 "Necropsy_Report_Printed" = col_character(),
                 "Data_transferred_to_MMLH_db" = col_character(),
                 "Photos_Uploaded" = col_character(),
                 "Photos_processed_renamed_deleted_if_not_useful" = col_character(),
                 "Photos_Printed" = col_character(),
                 "SPID_and_age_class_confirmed" = col_character(),
                 "Level_A_matches_SP_ID_and_age_class" = col_character(),
                 "SWFSC_RESPONSE" = col_character(),
                 "PickUp_Team_Contact" = col_character(),
                 "COMMENT" = col_character(),
                 "Biotox_sent" = col_character(),
                 "LAB_ID_for_genetics" = col_character())

# dat.cols <- c("text", "text", "text", "text", "text", "text", "text", "text", "text", "text",
#               "text", "text", "text", "text", "text", "text", "text", "text", "text", "text")

dat <- read_csv("data/DA EVENT DOLPHIN STRANDING LOG.csv",
                 col_types = dat.cols)

#.name_repair = function(col){ gsub(" ", "_", col) },

dat %>% 
  filter(SWFSC_RESPONSE == "Y") -> SW.resposne 

SW.resposne %>%
  