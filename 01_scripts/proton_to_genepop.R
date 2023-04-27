# Load and convert proton data to genepop format
# Sutherland Bioinformatics, 2022-09-09

proton_to_genepop <- function(hotspot_only = TRUE, neg_control="BLANK"){

  inputs <- list.files(path = "02_input_data/", pattern = ".xls")
  
  # Set nulls
  input.df <- NULL
  
  for(i in 1:length(inputs)){
    
    
    #### 01. Read in input, create run.name__barcode__sample.name identifier, then reduce to essential columns ####
    # Read in input
    print(paste0("Reading in file ", inputs[i]))
    input.df <- read.delim(file = paste0("02_input_data/", inputs[i]), header = T, sep = "\t")
    print("Data loaded; rows, cols:")
    print(dim(input.df))
    
    # Reduce columns
    print("Reducing columns to retain essentials")
    
    input.df  <- input.df[,c("Chrom", "Position", "Ref"
                              , "Variant", "Allele.Call", "Type", "Allele.Source"
                              , "Allele.Name", "Region.Name", "Coverage", "Strand.Bias"
                              , "Sample.Name", "Barcode", "Run.Name")]
    
    print("Columns limited; rows, cols:")
    print(dim(input.df))
    
    # Reporting
    print(paste0("Currently, there are ", length(unique(input.df$Allele.Name)), " unique markers"))
    num.markers <- length(unique(input.df$Allele.Name))
    
    
    # Remove non-hotspot if required
    if(hotspot_only==TRUE){
      
      # Reporting
      print("Removing novel (non-hotspot) markers")
      
      # Retain hotspot SNPs only
      input.df <- input.df[input.df$Allele.Source=="Hotspot",]
      
      print(paste0("Currently, there are ", length(unique(input.df$Allele.Name)), " unique markers"))
      
    }
    
    # Create new identifier comprised of Run Name, Barcode, and Sample Name
    input.df$identifier <- paste0(input.df$Run.Name, "__", input.df$Barcode, "__", input.df$Sample.Name)
    
    # Format into a matrix, genetic section
    input.df <- input.df[,c("identifier", "Allele.Name", "Ref", "Variant", "Allele.Call")]
    
    # Reporting
    print("Currently the dataset is comprised of: ")
    print(paste0(length(unique(input.df$identifier)), " samples"))
    print(paste0(length(unique(input.df$Allele.Name)), " markers"))
    
    # The data is now in 5 cols, with sample identifier, allele name, ref, variant, and the geno call for the sample
    head(input.df)
    
    #### 02. Convert allele calls to markers  ####
    # Per indiv, per marker, convert to actual genotype in genepop format
    
    # Create blank cols to be filled
    input.df$allele1 <- NA
    input.df$allele2 <- NA
    input.df$genepop <- NA
    
    # Reporting
    print("Converting Allele.Call to genotypes, incl. genepop format. Please be patient, this may take a while")
    print(paste0("Your dataset has ", nrow(input.df), " rows (1 row is 1 marker for 1 indiv.)"))
    
    # Loop to convert Allele.Call to allele 1 and 2 for each sample and marker
    print("Note: both actual nucleotides and the genepop format genotype will be retained in output")
    for(i in 1:nrow(input.df)){
      
      # Absent means homozygous reference
      if(input.df$Allele.Call[i]=="Absent"){
        
        input.df[i, "allele1"] <- input.df[i, "Ref"]
        input.df[i, "allele2"] <- input.df[i, "Ref"]
        input.df[i, "genepop"] <- "0101"
        
        # No Call means missing data  
      }else if(input.df$Allele.Call[i]=="No Call"){
        
        input.df[i, "allele1"] <- "missing"
        input.df[i, "allele2"] <- "missing"
        input.df[i, "genepop"] <- "0000"
        
        # Heterozygous means het 
      }else if(input.df$Allele.Call[i]=="Heterozygous"){
        
        input.df[i, "allele1"] <- input.df[i, "Ref"]
        input.df[i, "allele2"] <- input.df[i, "Variant"]
        input.df[i, "genepop"] <- "0102"
        
        # Homozygous means homozygous variant 
      }else if(input.df$Allele.Call[i]=="Homozygous"){
        
        input.df[i, "allele1"] <- input.df[i, "Variant"]
        input.df[i, "allele2"] <- input.df[i, "Variant"]
        input.df[i, "genepop"] <- "0202"
        
      }
    }
    
    # Reporting
    if(length(table(is.na(input.df$genepop)))==1){
      
      print("Completed conversion successfully")
      
    }else{
      
      #stop("There are NAs in the genepop formatted column, stopping function")
      ## the problem with this is that it will become hard to troubleshoot if it just crashes after the long run, #TODO: save output then crash
      
    }
    
    # Now have a complete, allele-based matrix
    head(input.df, n = 10)
    
    # # Save output
    # write.table(x = proton_trim.df, file = paste0("03_results/", proton.FN, "_proton_data_converted.txt")
    #             , quote = F, sep = "\t", row.names = F
    # )
     
    
    # Removing negative control samples
    # What are the neg. control samples?
    print(paste0("Removing, ", length(unique(input.df[grep(pattern = neg_control, x = input.df$identifier), "identifier"]))
                 , " negative control samples."))
    
    # Remove any neg. control samples
    print(paste0("Size of df prior to removal: ", nrow(input.df), " rows and ", ncol(input.df), " cols"))
    print(paste0("Includes ", length(unique(input.df$identifier)), " samples"))
    input.df <- input.df[grep(pattern = neg_control, x = input.df$identifier, invert = T), ]
    print(paste0("Size of df after removal: ", nrow(input.df), " rows and ", ncol(input.df), " cols"))
    print(paste0("Includes ", length(unique(input.df$identifier)), " retained samples"))
    
    ## Reporting
    # Output genotype summary
    print(table(input.df$genepop))
    print(paste0("Total number of genotypes in formatted output: ", sum(table(input.df$genepop))))
    
    
    
    
    
  }
  
}





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
