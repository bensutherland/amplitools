# Menu initiator for amplitools
# Input data should be tab-delimited text file with extension .xls from 
## Torrent VariantCaller output
# B. Sutherland, 2023-04-26

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Install packages


## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

file_sources <- list.files(path = "01_scripts/utilities/", pattern = "\\.r$", full.names = TRUE, ignore.case = TRUE)

# Source functions
for(fun in file_sources){
  print(fun)
  source(fun)
}
rm(fun, file_sources) # clean up

# Ready for analysis
options(scipen = 99999999)
