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
proton.FN <- "R_2022_08_04_09_19_56_user_S5XL-00533-1089-OYR-20220729_7 (2).xls"
hotspot_only <- TRUE # Set as true if want to only keep known SNPs

#### 01. Load and format data ####
input.df <- read.delim(file = paste0("02_input_data/", proton.FN), header = T, sep = "\t")
dim(input.df)
input.df[1:5,1:5]
colnames(input.df)

# Reduce to only keep necessary cols
proton.df  <- input.df[,c("Chrom", "Position", "Ref"
                      , "Variant", "Allele.Call", "Type", "Allele.Source"
                      , "Allele.Name", "Region.Name", "Coverage", "Strand.Bias"
                      , "Sample.Name", "Barcode", "Run.Name")]

rm(input.df) # Clean space

dim(proton.df)
head(proton.df)

# Backup
#proton.df.bck <- proton.df 

# Reporting
print(paste0("Currently, there are ", length(unique(proton.df$Allele.Name)), " unique markers"))

# Remove non-hotspot if required
if(hotspot_only==TRUE){
  
  # Reporting
  print("Removing novel (non-hotspot) markers")
  
  # Retain hotspot SNPs only
  proton.df <- proton.df[proton.df$Allele.Source=="Hotspot",]
  
}

#dim(proton.df)

print(paste0("Currently, there are ", length(unique(proton.df$Allele.Name)), " unique markers"))

# Create new identifier comprised of Run Name, Barcode, and Sample Name
proton.df$identifier <- paste0(proton.df$Run.Name, "__", proton.df$Barcode, "__", proton.df$Sample.Name)

# Format into a matrix, genetic section
proton_trim.df <- proton.df[,c("identifier", "Allele.Name", "Ref", "Variant", "Allele.Call")]
dim(proton_trim.df)
head(proton_trim.df)


#### 02. Convert from Allele.Call to actual markers ####
# Per indiv, per marker, provide true nucleotide genotype
#  Create new columns to be filled below
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
head(proton_trim.df)

# Save output
write.table(x = proton_trim.df, file = "proton_data_converted.txt", quote = F, sep = "\t", row.names = F)
#proton_trim.df.bck <- proton_trim.df


#### 03. Convert to genepop marker formats ####
proton_trim.df$genepop <- NA
head(proton_trim.df)


for(i in 1:nrow(proton_trim.df)){
  
  # Absent means homozygous reference
  if(proton_trim.df[i, "Allele.Call"]=="Absent"){
    
    proton_trim.df[i, "genepop"] <- "0101"
    
  # No Call means missing data
  }else if(proton_trim.df[i, "Allele.Call"]=="No Call"){
    
    proton_trim.df[i, "genepop"] <- "0000"
    
  # Heterozygous means heterozygous
  }else if(proton_trim.df[i, "Allele.Call"]=="Heterozygous"){
    
    proton_trim.df[i, "genepop"] <- "0102"
    
  # Homozygous means homozygous alternate
  }else if(proton_trim.df[i, "Allele.Call"]=="Homozygous"){
    
    proton_trim.df[i, "genepop"] <- "0202"
    
  }
}

# See summary of data
table(proton_trim.df$genepop)
#proton_trim.df.bck2 <- proton_trim.df

#### 03. Structural change to matrix to prepare for genepop format ####

# REMOVE THE BLANKs
unique(proton_trim.df$Sample.Name)
dim(proton_trim.df)
head(proton_trim.df)

# backup
#proton.df.bck.bck <- proton_trim.df

proton_trim.df <- proton_trim.df[proton_trim.df$Sample.Name!="BLANK", ]
dim(proton_trim.df)

head(proton_trim.df)
str(proton_trim.df)

# Identify all samples present
samples <- unique(proton_trim.df$Sample.Name)

# Do all of the samples have an equal number of records? 
table(table(proton_trim.df$Sample.Name)==592)
### TODO: This can be improved ###

# Debugging
#samples <- head(samples, n = 110)

# Set nulls
soi <- NULL; line.item <- NULL; genetics.df <- NULL; line.item.df <- NULL

for(i in 1:length(samples)){
  
  print(paste0("i = ", i))
  
  # Take each sample
  soi <- samples[i]
  print(soi)
  
  # Take all of the genepop identifiers
  line.item <- proton_trim.df[proton_trim.df$Sample.Name==soi, c("Allele.Name", "genepop")]
  #head(line.item)
  
  # Set colnames using Allele.Name, in horizontal form now? 
  colnames(x = line.item)[which(colnames(line.item)=="genepop")] <- soi
  
  # Make a df
  line.item.df <- as.data.frame(line.item)
  
  #head(line.item.df)
  
  # Build a full df
  if(i==1){
    
    genetics.df <- line.item.df
    
  }else if(i>1){
    
    genetics.df <- merge(x = genetics.df, y = line.item.df, by = "Allele.Name", all.x = T)
    # SHOULD INCLUDE ALL REQUIREMENT
    
  }
  
  print(head(genetics.df))
  
}

dim(genetics.df)

genetics.df[1:10, 1:10]



write.table(x = genetics.df, file = "genetic_data_only.txt", quote = F, sep = "\t", row.names = F)

genetics_prep.df <- t(genetics.df)
genetics_prep.df[1:10, 1:10]
dim(genetics_prep.df)
indiv <- rownames(genetics_prep.df)
test <- cbind(rownames(genetics_prep.df), genetics_prep.df)
test[1:5,1:5]

colnames(x = test) <- test[1,]
test[1:5,1:5]

colnames(test)[which(colnames(test)=="Allele.Name")] <- "indiv"
test[1:5,1:5]

# Keep all others but not that row
test2 <- test[grep(pattern = "Allele.Name", x = test[,"indiv"], invert = T), ]

test2[1:10,1:10]

write.table(x = test2, file = "genetic_data_only_final.txt", quote = F, sep = "\t", row.names = F)

# Move to shell to finalize the genepop

