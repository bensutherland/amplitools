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
neg_control <- "BLANK" # Set based on the name of the negative control wells in the study

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
num.markers <- length(unique(proton.df$Allele.Name))

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


#### 02. Convert from Allele.Call to actual markers, incl. genepop format ####
# Per indiv, per marker, provide true nucleotide genotype (and genepop format)
#  Create new columns to be filled below
proton_trim.df$allele1 <- NA
proton_trim.df$allele2 <- NA
proton_trim.df$genepop <- NA

# # For debugging
# proton_trim.df <- head(proton_trim.df, n = 200)
# head(proton_trim.df)

print("Converting Allele.Call to genotypes, incl. genepop format. Please be patient, this may take a while")

# Loop to convert Allele.Call to allele 1 and 2 for each sample and marker
for(i in 1:nrow(proton_trim.df)){
  
  # Absent means homozygous reference
  if(proton_trim.df$Allele.Call[i]=="Absent"){
    
    proton_trim.df[i, "allele1"] <- proton_trim.df[i, "Ref"]
    proton_trim.df[i, "allele2"] <- proton_trim.df[i, "Ref"]
    proton_trim.df[i, "genepop"] <- "0101"
    
  # No Call means missing data  
  }else if(proton_trim.df$Allele.Call[i]=="No Call"){
    
    proton_trim.df[i, "allele1"] <- "missing"
    proton_trim.df[i, "allele2"] <- "missing"
    proton_trim.df[i, "genepop"] <- "0000"
  
  # Heterozygous means het 
  }else if(proton_trim.df$Allele.Call[i]=="Heterozygous"){
    
    proton_trim.df[i, "allele1"] <- proton_trim.df[i, "Ref"]
    proton_trim.df[i, "allele2"] <- proton_trim.df[i, "Variant"]
    proton_trim.df[i, "genepop"] <- "0102"
    
  # Homozygous means homozygous variant 
  }else if(proton_trim.df$Allele.Call[i]=="Homozygous"){
    
    proton_trim.df[i, "allele1"] <- proton_trim.df[i, "Variant"]
    proton_trim.df[i, "allele2"] <- proton_trim.df[i, "Variant"]
    proton_trim.df[i, "genepop"] <- "0202"
    
  }
}

# Now have a complete, allele-based matrix
head(proton_trim.df, n = 10)

# Save output
write.table(x = proton_trim.df, file = paste0("03_results/", proton.FN, "_proton_data_converted.txt")
            , quote = F, sep = "\t", row.names = F
            )

# Backup
#proton_trim.df.bck <- proton_trim.df

# What are the neg. control samples? 
unique(proton_trim.df[grep(pattern = neg_control, x = proton_trim.df$identifier), "identifier"])

# Remove any neg. control samples
dim(proton_trim.df)
proton_trim.df <- proton_trim.df[grep(pattern = neg_control, x = proton_trim.df$identifier, invert = T), ]
dim(proton_trim.df)
head(proton_trim.df)

# See summary of data
table(proton_trim.df$genepop)


#### 03. Change structure into matrix to prepare for genepop format ####
# Confirm that all of the samples have an equal number of records
table(table(proton_trim.df$identifier)==num.markers)

# Identify all samples present
samples <- unique(proton_trim.df$identifier)

# Per sample, loop to collect the marker names and genepop genotypes, now labeled with identifier
# Set nulls
soi <- NULL; line.item <- NULL; genetics.df <- NULL; line.item.df <- NULL

for(i in 1:length(samples)){
  
  # Identify the sample name for this round
  soi <- samples[i]
  
  # Reporting
  print(paste0("i = ", i))
  print(soi)
  
  # For this sample, retain all marker names and corresponding genepop genotypes
  line.item <- proton_trim.df[proton_trim.df$identifier==soi, c("Allele.Name", "genepop")]
  
  # For this record, set the genepop column name as the sample name
  colnames(x = line.item)[which(colnames(line.item)=="genepop")] <- soi
  
  # Convert record to df
  line.item.df <- as.data.frame(line.item)
  
  # Build a full df using the record
  # If it is the first record, initialize the df
  if(i==1){
    
    genetics.df <- line.item.df
    
  # If it is a subsequent record, merge with the existing df by marker names
  }else if(i>1){
    
    genetics.df <- merge(x = genetics.df, y = line.item.df, by = "Allele.Name", all.x = T)
    # Note: the all.x is necessary here to avoid issues if the number of markers is not identical
    
  }
  
  #print(head(genetics.df))
  
}

# Observe output
dim(genetics.df)
genetics.df[1:4, 1:4]

# Save output
write.table(x = genetics.df, file = paste0("03_results/", proton.FN, "_genetic_data_only.txt"), quote = F, sep = "\t", row.names = F)

# Reshape the output; transpose to columns as markers
genetics_prep.df <- t(genetics.df)
genetics_prep.df[1:4, 1:4]
dim(genetics_prep.df)


# Include the rownames (indiv names) as a new column
genotypes.df <- cbind(rownames(genetics_prep.df), genetics_prep.df)
genotypes.df[1:3,1:3]

# Use the marker names (the first row) as colnames for the object
colnames(x = genotypes.df) <- genotypes.df[1,]
genotypes.df[1:5,1:5]

# Rename the column 'Allele.Name' as the more correct, 'indiv' 
colnames(genotypes.df)[which(colnames(genotypes.df)=="Allele.Name")] <- "indiv"
genotypes.df[1:3,1:3]

# Remove the row containing the old header information (used to contain markers)
genotypes.df <- genotypes.df[grep(pattern = "Allele.Name", x = genotypes.df[,"indiv"], invert = T), ]

# View
genotypes.df[1:5,1:5]

# Save output
paste0("03_results/", proton.FN, "_genetic_data_only.txt")
write.table(x = genotypes.df, file = paste0("03_results/", proton.FN, "_genetic_data_only_final.txt"), quote = F, sep = "\t", row.names = F)

# See README for directions on completing the genepop in shell
