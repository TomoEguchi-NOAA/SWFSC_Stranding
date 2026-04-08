# SWFSC_Stranding_fcns.R
# 
# Functions common to many scripts in the SWFSC_Stranding project

connection.string <- function(database){
  return(paste0("Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-s;Database=",
                database, ";Trusted_Connection=yes;TrustServerCertificate=yes;"))
  
  # return(paste0("Driver={ODBC Driver 18 for SQL Server};Server=swc-estrella-ut.nmfs.local;Database=",
  #               database, ";Trusted_Connection=yes;Port=1433;TrustServerCertificate=yes;"))

  }

