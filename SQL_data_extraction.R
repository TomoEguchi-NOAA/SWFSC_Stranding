# Access stranding and life history databases and extract morphological and
# life history data. This project is for a Hollings scholar (Raija Hammond)
# in Summer 2026.
# 
rm(list = ls())
library(RODBC)
library(tidyverse)
library(readr)

MMLH.con <- "Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-s;Database=MMLH;Trusted_Connection=yes;TrustServerCertificate=yes;"
MMLH.2019 <- odbcDriverConnect(MMLH.con)
MMLH.tables <- sqlTables(MMLH.2019)

Common.con <- "Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-s;Database=SWFSCCommon;Trusted_Connection=yes;TrustServerCertificate=yes;"
Common  <- odbcDriverConnect(Common.con)

Tissue.con <- "Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-s;Database=TissueArchive;Trusted_Connection=yes;TrustServerCertificate=yes;"
Tissue  <- odbcDriverConnect(Tissue.con)

odbcCloseAll()
