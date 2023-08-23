# Load and convert proton data to genepop format
# Sutherland Bioinformatics, 2022-09-09

proton_to_genepop <- function(hotspot_only = TRUE, neg_control="BLANK"){

  # Load vc data into input.list
  load_vc(input_folder = "02_input_data")
  
  input.df <- NULL
  for(c in 1:length(names(input.list))){
    
    input.df <- input.list[[c]]
    
    # Remove non-hotspot if required
    if(hotspot_only==TRUE){
      
      # Reporting
      print("Removing novel (non-hotspot) markers")
      
      # Retain hotspot SNPs only
      input.df <- input.df[input.df$Allele.Source=="Hotspot",]
      
      print(paste0("Currently, there are ", length(unique(input.df$Allele.Name)), " unique markers"))
      
    # If hotspot only is not true, then keep all SNP variants
    }else if(hotspot_only==FALSE){
      
      # Reporting
      print("Keeping all SNP variants, including non-hotspot variants")
      
      # Stop execution, this function has not yet been developed
      stop("**Process stopped: characterizing non-hotspot variants (i.e., 'novel SNPs') is not yet implemented, but is planned for future development**")
      
      # Retain all SNPs
      #input.df <- input.df[input.df$Type=="SNP", ]
      
      #print(paste0("Currently, there are ", length(unique(input.df$Allele.Name)), " unique markers"))
      
      #if(keep_type=="top_maf"){
      #  }
      
    }
    

    ## Summarize data characteristics
    # Per sample mean marker read depth
    print("Summarizing per sample mean/median marker read depth")
    mean_cov_per_sample.df   <- aggregate(Coverage ~ identifier, data = input.df, mean)
    median_cov_per_sample.df <- aggregate(Coverage ~ identifier, data = input.df, median)
    
    coverage_fig.FN <- paste0("03_results/", "per_sample_marker_depth_", gsub(pattern = "\\.df", x = names(input.list)[c], replacement = ""), ".pdf")
    
    pdf(file = coverage_fig.FN, width = 5, height = 7)
    par(mfrow=c(2,1))
    hist(x = mean_cov_per_sample.df$Coverage, breaks = 20, main = ""
         , xlab = "Per sample mean marker read depth"
    )
    
    hist(x = median_cov_per_sample.df$Coverage, breaks = 20, main = ""
         , xlab = "Per sample median marker read depth"
    )
    dev.off()
    rm(coverage_fig.FN)
    
    # Summarize per marker mean marker read depth (summarize across samples)
    print("Summarizing per marker mean/median read depth, summarizing across samples")
    mean_cov_per_marker.df   <- aggregate(Coverage ~ Allele.Name, data = input.df, mean)
    median_cov_per_marker.df <- aggregate(Coverage ~ Allele.Name, data = input.df, median)
    
    coverage_fig.FN <- paste0("03_results/", "per_marker_depth_across_samples_", gsub(pattern = "\\.df", x = names(input.list)[c], replacement = ""), ".pdf")
    pdf(file = coverage_fig.FN, width = 5, height = 7)
    par(mfrow=c(2,1))
    hist(x = mean_cov_per_marker.df$Coverage, breaks = 20, main = ""
         , xlab = "Marker average read depth (summarized across samples)"
    )
    
    hist(x = median_cov_per_marker.df$Coverage, breaks = 20, main = ""
         , xlab = "Marker median read depth (summarized across samples)"
    )
    dev.off()
    
    
    # Reporting
    print(paste0(
      "Negative control samples have an average mean marker coverage of "
      , round(mean(mean_cov_per_sample.df[grep(pattern = neg_control, x = mean_cov_per_sample.df$identifier), "Coverage"]), digits = 3)
    ))
    
    print(paste0(
      "Experimental samples have an average mean marker coverage of "
      , round(mean(mean_cov_per_sample.df[grep(pattern = neg_control, x = mean_cov_per_sample.df$identifier, invert = T), "Coverage"]), digits = 3)
    ))
    
    # Save out data characteristics
    coverage.FN <- paste0("03_results/", "mean_marker_cov_per_sample_", gsub(pattern = "\\.df", x = names(input.list)[c], replacement = ""), ".txt")
    write.table(x = mean_cov_per_sample.df, file = coverage.FN, quote = F, col.names = T, row.names = F, sep = "\t")
    
    # Save out data characteristics
    coverage.FN <- paste0("03_results/", "mean_marker_cov_per_marker_", gsub(pattern = "\\.df", x = names(input.list)[c], replacement = ""), ".txt")
    write.table(x = mean_cov_per_marker.df, file = coverage.FN, quote = F, col.names = T, row.names = F, sep = "\t")
    
    
    ### Prepare data into matrix
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
    
    ### TODO: output this as a data output (nucleotides lost below)
    
    # # Save output
    # write.table(x = proton_trim.df, file = paste0("03_results/", proton.FN, "_proton_data_converted.txt")
    #             , quote = F, sep = "\t", row.names = F
    # )
    
    
    #### 03. Restructure to create genotype matrix for genepop format  ####
    ## Data check
    if(length(table(table(input.df$identifier)))==1){
      
      print("All individuals have an equal number of rows (markers), as required, in total per indiv: ")
      print(length(unique(input.df$Allele.Name)))
      
    }else{
      
      print("Warning! Not all individuals have the equal number of rows, and so something has gone wrong. Need to debug.")
      
    }
      
    
    # Reporting
    print("Building genotype matrix (df) to be used for genepop construction")
    
    # Identify all samples present
    samples <- unique(input.df$identifier)
    
    ## Per sample, collect the marker names and genepop-formatted genotypes
    # Set nulls
    soi <- NULL; line.item <- NULL; genetics.df <- NULL; line.item.df <- NULL
    # Loop
    for(j in 1:length(samples)){
      
      # Identify the sample name for this round
      soi <- samples[j]
      
      ## Reporting (debugging)
      # print(paste0("j = ", j))
      # print(soi)
      
      # For this sample, retain all marker names and corresponding genepop genotypes
      line.item <- input.df[input.df$identifier==soi, c("Allele.Name", "genepop")]
      
      # Set the genepop genotype column name as the sample name
      colnames(x = line.item)[which(colnames(line.item)=="genepop")] <- soi
      
      # Convert to df
      line.item.df <- as.data.frame(line.item)
      rm(line.item)
      
      # Build the full genotype df
      # If this is the first record, initialize the df
      if(j==1){
        
        genetics.df <- line.item.df
        
        # If it is a subsequent record, merge with the existing df by marker names
      }else if(j>1){
        
        genetics.df <- merge(x = genetics.df, y = line.item.df, by = "Allele.Name", all.x = T)
        # Note: the all.x is necessary here to avoid issues if the number of markers is not identical
        
      }
      
      #print(head(genetics.df))
      
    }
    
    # Observe output
    dim(genetics.df)
    genetics.df[1:4, 1:4]
    
    # # Save output
    # write.table(x = genetics.df, file = paste0("03_results/", proton.FN, "_genetic_data_only.txt"), quote = F, sep = "\t", row.names = F)
    
    # Reshape the output; transpose to make samples as rows and markers as columns
    genetics_prep.df <- t(genetics.df)
    genetics_prep.df[1:4, 1:4]
    print(dim(genetics_prep.df))
    rm(genetics.df)
    
    # Include the rownames (indiv names) as a new column
    genotypes.df <- cbind(rownames(genetics_prep.df), genetics_prep.df)
    genotypes.df[1:3,1:3]
    rm(genetics_prep.df)
    
    # Use the marker names (the first row) as colnames for the object
    colnames(x = genotypes.df) <- genotypes.df[1,]
    genotypes.df[1:5,1:5]
    
    # Rename the column 'Allele.Name' as the more correct, 'indiv' 
    colnames(genotypes.df)[which(colnames(genotypes.df)=="Allele.Name")] <- "indiv"
    genotypes.df[1:3,1:3]
    
    # Remove the row containing the old header information
    genotypes.df <- genotypes.df[grep(pattern = "Allele.Name", x = genotypes.df[,"indiv"], invert = T), ]
    
    # View
    genotypes.df[1:5,1:5]
    
    # Reporting
    print("The genotype matrix is now prepared.")
    
    # Create output filename
    output.FN <- paste0("02_input_data/prepped_matrices/"
                        , gsub(pattern = "\\.df", replacement = "", x = names(input.list)[c])
                        , "_gen_data.txt"
                        )
    
    # Reporting
    print(paste0("Saving output to ", output.FN))
    
    # Save output
    write.table(x = genotypes.df, file = output.FN, quote = F, sep = "\t", row.names = F)
    
  }
  
}

# See README for directions on completing the genepop in shell
