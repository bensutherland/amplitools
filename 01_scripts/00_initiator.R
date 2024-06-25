# Menu initiator for amplitools
# Input data should be tab-delimited text file with extension .xls from 
## Torrent VariantCaller output
# B. Sutherland, 2023-04-26

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Install and load packages
# Devtools
#  Mac instructions:  
#   #install.packages("devtools")
#   #devtools::install_github("hadley/devtools")
# Ubuntu suggest installing devtools with apt-get or similar package manager

#install.packages("Rcpp")
#install.packages("tidyverse")
#devtools::install_github("eriqande/CKMRsim", build_vignettes = TRUE)
#install.packages("igraph")
#install.packages("adegenet")
#devtools::install_github("delomast/EFGLmh") # required for microhaplotype workflow

#library("Rcpp")
library("devtools")
library("tidyverse")
library("CKMRsim")
library("igraph")
library("adegenet")
library("EFGLmh") # required for microhaplotype workflow

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

print("amplitools is loaded and ready for analysis")
