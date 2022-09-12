# Load and convert proton data to genepop format
# Sutherland Bioinformatics, 2022-09-09

# clear workspace
#rm(list=ls())

#### Front Matter ####
options(scipen = 99999999)

## Install packages and load libraries


## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

## User set variables
proton.FN <- "~/Documents/00_sutherland_bioinformatics/GBMF_UBC_Pacific_oyster/amplicon_panel/00_results_of_pilot_study_2022-09-07/R_2022_08_04_09_19_56_user_S5XL-00533-1089-OYR-20220729_7 (2).xls"
hotspot_only <- TRUE # Set as true if want to only keep known SNPs

#### 00. Load data ####
input.df <- read.delim(file = proton.FN, header = T, sep = "\t")
dim(input.df)
input.df[1:5,1:5]
colnames(input.df)

# Reduce to only keep necessary cols
proton.df  <- input.df[,c("Chrom", "Position", "Ref"
                      , "Variant", "Allele.Call", "Type", "Allele.Source"
                      , "Allele.Name", "Region.Name", "Coverage", "Strand.Bias"
                      , "Sample.Name", "Barcode", "Run.Name")]

rm(input.df) # Clean space

head(proton.df)

# Backup
#proton.df.bck <- proton.df 

dim(proton.df)

print(paste0("There are ", length(unique(proton.df$Allele.Name)), " unique markers"))

# Remove non-hotspot if required
if(hotspot_only==TRUE){
  
  # Reporting
  print("Removing novel (non-hotspot) markers")
  
  # Retain hotspot SNPs only
  proton.df <- proton.df[proton.df$Allele.Source=="Hotspot",]
  
}

#dim(proton.df)

print(paste0("There are ", length(unique(proton.df$Allele.Name)), " unique markers"))

# Format into a matrix, genetic section
proton_trim.df <- proton.df[,c("Sample.Name", "Allele.Name", "Ref", "Variant", "Allele.Call")]
head(proton_trim.df)

#TODO: Could add additional QC here

#### 01. Convert from Allele.Call to actual markers ####
# Assign per indiv, per marker true markers for each allele
# Create new columns
proton_trim.df$allele1 <- NA
proton_trim.df$allele2 <- NA

# # For debugging
# proton_trim.df <- head(proton_trim.df, n = 200)
# head(proton_trim.df)

# Loop to convert Allele.Call to allele 1 and 2 for each sample and marker
for(i in 1:nrow(proton_trim.df)){
  
  # Absent means homozygous reference
  if(proton_trim.df$Allele.Call[i]=="Absent"){
    
    proton_trim.df[i, "allele1"] <- proton_trim.df[i, "Ref"]
    proton_trim.df[i, "allele2"] <- proton_trim.df[i, "Ref"]
  
  # No Call means missing data  
  }else if(proton_trim.df$Allele.Call[i]=="No Call"){
    
    proton_trim.df[i, "allele1"] <- "missing"
    proton_trim.df[i, "allele2"] <- "missing"
  
  # Heterozygous means het 
  }else if(proton_trim.df$Allele.Call[i]=="Heterozygous"){
    
    proton_trim.df[i, "allele1"] <- proton_trim.df[i, "Ref"]
    proton_trim.df[i, "allele2"] <- proton_trim.df[i, "Variant"]
  
  # Homozygous means homozygous variant 
  }else if(proton_trim.df$Allele.Call[i]=="Homozygous"){
    
    proton_trim.df[i, "allele1"] <- proton_trim.df[i, "Variant"]
    proton_trim.df[i, "allele2"] <- proton_trim.df[i, "Variant"]
    
  }
}

# Now have a complete, allele-based matrix



















