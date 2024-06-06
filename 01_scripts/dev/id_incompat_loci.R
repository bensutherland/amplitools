# Identify Mendelian incompatibility loci based on empirical data
# 2024-06-07; Ben J. G. Sutherland (VIU, SBIO)
# Requires that 01_scripts/utilities/ckmr_from_rubias.R has been run

id_incompat_loci <- function(data = obj, input_report.FN = "03_results/ckmr_input_rubias_mhaps_396_2024-05-24_VIU_F2_vs_VIU_F1_2024-06-05/po_VIU_F1_vs_VIU_F2_pw_logl_5_report.txt"
                             , strong_logl_cutoff = 10, strong_parent_obs = 2, num_offsp_incompat = 4){
 
  #### Identify strong empirical trios as 'true' relationships ####
  ## Load ckmrsim report produced by ckmr_from_rubias()
  po_report.df <- read.table(file = input_report.FN, header = T, sep = "\t")
  head(po_report.df)
  dim(po_report.df)
  po_report.df$p1_logl <- as.numeric(po_report.df$p1_logl) # ensure numeric
  po_report.df$p2_logl <- as.numeric(po_report.df$p2_logl) # ensure numeric
  
  ## Limit report df to only "strong IDs"
  # Remove records where one parent is missing
  po_report.df <- po_report.df[!is.na(po_report.df$p1_logl) & !is.na(po_report.df$p2_logl), ] 
  dim(po_report.df)
  head(po_report.df)
  
  # Keep trios where P-O logls are greater than cutoff, and remove any trios with more than two parents assigned above the cutoff
  po_report.df <- po_report.df[po_report.df$p1_logl >= strong_logl_cutoff & po_report.df$p2_logl >= strong_logl_cutoff & is.na(x = po_report.df$other_assigns), ]
  dim(po_report.df)
  head(po_report.df)
  
  ## Create empirical family map
  # Identify all parental pairs from the strong trios, ordered alphabetically within a parental pair
  parents <- NULL; parent.vec <- NULL
  for(i in 1:nrow(po_report.df)){
    
    parents <- sort(c(po_report.df$p1[i], po_report.df$p2[i]))
    
    parents <- paste0(parents, collapse = "__")
    
    # Build the parent vector
    parent.vec <- c(parent.vec, parents)
    
  }
  rm(parents) # clean up
  
  empirical_family_map.df <- as.data.frame(parent.vec)
  rm(parent.vec)
  colnames(empirical_family_map.df)[1] <- "family.id"
  head(empirical_family_map.df)
  
  # Add other cols (note: requires no sort has occurred)
  empirical_family_map.df$offspring <- po_report.df$indiv
  empirical_family_map.df <- separate(data = empirical_family_map.df, col = "family.id", into = c("parent1", "parent2"), sep = "__", remove = F)
  
  empirical_family_map.df <- empirical_family_map.df[,c("offspring", "family.id", "parent1", "parent2")]
  head(empirical_family_map.df)
  
  # Determine the frequency of parental pair observations
  family_freq.df <- as.data.frame(table(empirical_family_map.df$family.id))
  dim(family_freq.df)
  
  # Only keep those with more than the cutoff number observation
  family_freq.df <- family_freq.df[family_freq.df$Freq >= strong_parent_obs, ]
  dim(family_freq.df)
  family_freq.df
  select_empirical_families <- family_freq.df$Var1
  
  # Select only those families from the family map, this will comprise the empirical trios
  empirical_family_map.df <- empirical_family_map.df[empirical_family_map.df$family.id %in% select_empirical_families, ]
  empirical_family_map.df
  dim(empirical_family_map.df)
  length(unique(empirical_family_map.df$family.id))
  unique(empirical_family_map.df$family.id)
  
  
  #### Back up prepared obj ####
  #obj.bck <- data
  
  
  #### 01.b. Create genind object with only the offspring and parents from the strong trios ####
  offspring_to_keep <- empirical_family_map.df$offspring
  parents_to_keep   <- unique(c(empirical_family_map.df$parent1,  empirical_family_map.df$parent2))
  inds_to_keep <- c(offspring_to_keep, parents_to_keep)
  length(offspring_to_keep)
  length(parents_to_keep)
  length(inds_to_keep)
  
  data_strong <- data[inds_to_keep]
  data_strong
  
  #### THIS COULD BE A STOP POINT FOR A FUNCTION ####
  
  
  #### 02. Determine expected offspring genos ####
  # What are the unique parent combinations? 
  head(empirical_family_map.df)
  dim(empirical_family_map.df)
  unique_families.df <- empirical_family_map.df[!duplicated(x = empirical_family_map.df$family.id), ]
  unique_families.df <- unique_families.df[, c("family.id", "parent1", "parent2")]
  head(unique_families.df)
  
  # Create empty matrix to fill with possible offspring genos based on parents
  loci <- matrix(data = NA, nrow = nrow(unique_families.df), ncol = nLoc(obj))
  colnames(loci) <- locNames(obj)
  loci[1:5,1:5]
  unique_fams_and_loci.df <- cbind(unique_families.df, loci)
  unique_fams_and_loci.df <- as.data.frame(unique_fams_and_loci.df)
  unique_fams_and_loci.df[1:5,1:5] # this will be filled
  dim(unique_fams_and_loci.df)
  
  #### TODO: this should be a separate function ####
  ## Infer the potential offspring genotypes per pairing
  # set nulls
  p1 <- NULL; p2 <- NULL; fam <- NULL; family_marker.list <- list()
  for(i in 1:nrow(unique_fams_and_loci.df)){
    
    # Set variables
    fam <- unique_fams_and_loci.df$family.id[i]
    fam <- as.character(fam) # necessary in case it is being interpreted as an integer
    p1 <-  unique_fams_and_loci.df$parent1[i]
    p2 <-  unique_fams_and_loci.df$parent2[i]
    
    # Reporting
    print(paste0("Analyzing empirical family ", fam, ", comprised of parents: ", p1, " and ", p2))
    
    # Subset the genind to only the parents for the family
    parents.obj <- data_strong[c(p1, p2)]
    
    loi <- NULL; p1.alleles <- NULL; p2.alleles <- NULL; offspring_genos <- NULL
    opt1 <- NULL; opt2 <- NULL; opt3 <- NULL; opt4 <- NULL
    for(l in 1:nLoc(data_strong)){
      
      # Set locus of interest
      loi <- locNames(data_strong)[l]
    
      # Confirm that neither parent is missing, otherwise indicate
      !isTRUE(all(is.na(data_strong[p1, loc = loi]$tab)))
      !isTRUE(all(is.na(data_strong[p2, loc = loi]$tab)))
      
      !isTRUE(all(is.na(data_strong[p1, loc = loi]$tab))) # not missing
      !isTRUE(all(is.na(data_strong[p2, loc = loi]$tab))) # is missing
      
      # if neither parent is NA for the locus, infer potential genotypes
      if(!is.na(data_strong[p1, loc = loi]$tab[1]) & !is.na(data_strong[p2, loc = loi]$tab[1])){
        
        # Obtain the parental alleles
        p1.alleles <- colnames(data_strong[p1, loc = loi]$tab)[which(data_strong[p1, loc = loi]$tab > 0)]
        p2.alleles <- colnames(data_strong[p2, loc = loi]$tab)[which(data_strong[p2, loc = loi]$tab > 0)]
        
        # if either parent only is showing a single allele, assume it is homozygous for this allele (i.e., diploid is assumed)
        if(length(p1.alleles)==1){
          p1.alleles <- c(p1.alleles, p1.alleles)
        }
        
        if(length(p2.alleles)==1){
          p2.alleles <- c(p2.alleles, p2.alleles)
        }
        
        # What are the putative offspring genotypes? Will need to sort the offspring too for matching purposes
        opt1 <- sort(c(p1.alleles[1], p2.alleles[1]))
        opt1 <- paste0(opt1, collapse = ";")
        
        opt2 <- sort(c(p1.alleles[1], p2.alleles[2]))
        opt2 <- paste0(opt2, collapse = ";")
        
        opt3 <- sort(c(p1.alleles[2], p2.alleles[1]))
        opt3 <- paste0(opt3, collapse = ";")
        
        opt4 <- sort(c(p1.alleles[2], p2.alleles[2]))
        opt4 <- paste0(opt4, collapse = ";")
        
        offspring_genos <- paste0(opt1, ";;", opt2, ";;", opt3, ";;", opt4)
        
        unique_fams_and_loci.df[unique_fams_and_loci.df$family.id==fam, loi] <- offspring_genos
        
        
      }else{
        
        # set this as missing data
        offspring_genos <- paste0(NA, ";;", NA, ";;", NA, ";;", NA)
        
        unique_fams_and_loci.df[unique_fams_and_loci.df$family.id==fam, loi] <- offspring_genos
        
        
      }
      
      
    }
  }
    
    unique_fams_and_loci.df[1:5, 1:5]
   
  
    # Save out parental genos per marker
    write.csv(x = unique_fams_and_loci.df, file = "03_results/parental_genos_per_marker.csv"
              , quote = F, row.names = F
              )
    
    

  #### 03. Evaluate offspring against parents ####
  #### TODO: this should be a separate function ####
    
  #### Needs a rethink #####
    
  # Who are the offspring of the family? 
  indNames(data_strong) # individuals available
  head(empirical_family_map.df)
  
  families <- unique(empirical_family_map.df$family.id)
  
  fam <- NULL; p1 <- NULL; p2 <- NULL; keep <- NULL; offspring <- NULL
  for(o in 1:length(families)){
    
    # Set family
    fam <- families[o]
    p1 <- gsub(pattern = "__.*", replacement = "", x = fam)
    p2 <- gsub(pattern = ".*__", replacement = "", x = fam)
    
    # Reporting
    print(paste0("Analyzing family ", fam, ", comprised of parents: ", p1, " and ", p2))
    print(paste0("Empirically defined offspring are selected from the empirical family map"))
    
    # Retain only the offspring empirically determined to be offspring of this specific family
    keep <- empirical_family_map.df[empirical_family_map.df$family.id==fam, "offspring"]
    offspring <- data_strong[keep]
    print(indNames(offspring))
    
    offspring_genos.mat <- matrix(data = NA, nrow = nrow(empirical_family_map.df), ncol = nLoc(offspring))
    colnames(offspring_genos.mat) <- locNames(offspring)
    empirical_family_map_w_scores.df <- cbind(empirical_family_map.df, offspring_genos.mat)
    empirical_family_map_w_scores.df[1:5,1:5]
    
    # Loop over to check the offspring's results
    ooi <- NULL
    for(i in 1:nInd(offspring)){
      
      # Define the offspring
      ooi <- indNames(offspring)[i]
      
      # Define the family
      foi <- empirical_family_map.df[empirical_family_map.df$offspring==ooi, "family.id"] # technically not necessary, we define above
      
      loi <- NULL; offspring_geno <- NULL; opts <- NULL
      for(l in 1:nLoc(offspring)){
        
        # Define the locus
        loi <- locNames(offspring)[l]
        
        # What is the offspring's geno for this locus?
        # THIS IS AN ISSUE IF HOMOZYGOUS ### TODO ####
        offspring_geno <- paste0(sort(colnames(offspring[ooi, loc=loi]$tab)[which(offspring[ooi, loc=loi]$tab > 0)]), collapse = ";")
        
        # Look-up in parental df
        opts <- strsplit(x = unique_fams_and_loci.df[unique_fams_and_loci.df$family.id==foi, loi], split = ";;")
        opts <- unlist(opts)
        
        # If this locus is missing data in the parents
        if(all((unlist(strsplit(x = opts, split = ";;"))=="NA"))){
          
          empirical_family_map_w_scores.df[empirical_family_map_w_scores.df$offspring==ooi, loi] <- "missing"
          
          # If the locus is not missing in the parents
        }else{
          
          # if it is present in parental opts, give it a TRUE
          if(offspring_geno %in% opts){
            
            empirical_family_map_w_scores.df[empirical_family_map_w_scores.df$offspring==ooi, loi] <- "OK"
            
          }else if(!offspring_geno %in% opts){
            
            empirical_family_map_w_scores.df[empirical_family_map_w_scores.df$offspring==ooi, loi] <- offspring_geno
            
          }
        }
      }
    }
  }
    
  
  empirical_family_map_w_scores.df[1:5,1:5]
  

  # Save out
  write.csv(x = empirical_family_map_w_scores.df, file = "03_results/offspring_compats_per_marker.csv"
            , quote = F, row.names = F
  )
  
  
  #### WORKING HERE ####
  
    
    
    
    
    
    
    
    
    
    #### WON'T WORK HERE ####
    # Create empty df
    single_genos.df <- NULL
    
    # For the offspring genotypes, remove the second allele column, keeping only the first allele column
    single_genos.df <- offspring$tab[, grep(pattern = "\\.02", x = colnames(x = offspring$tab), invert = T)]
    #colnames(single_genos.df)
    
    # For these offspring, edit the single_genos df to give the genotype as calculated from the single col data
    for(n in 1:nrow(single_genos.df)){
      
      single_genos.df[n,] <- gsub(pattern = "1", replacement = "het", x = single_genos.df[n,])
      single_genos.df[n,] <- gsub(pattern = "2", replacement = "homo.ref", x = single_genos.df[n,])
      single_genos.df[n,] <- gsub(pattern = "0", replacement = "homo.alt", x = single_genos.df[n,])
      
    }
    
    #single_genos.df[,1:10]
    
    # Write out offspring geno calls
    write.csv(x = single_genos.df, file = paste0("03_results/offspring_geno_calls_", fam, ".csv"))
    
    ## Check the observed offspring genotypes against the expected offspring genotypes, derived from the parents
    # Prepare a T/F matrix to fill
    single_genos_eval.df <- single_genos.df
    
    mname <- NULL
    for(i in 1:ncol(single_genos.df)){
      
      # Check each marker
      mname <- colnames(single_genos.df)[i]
      
      # What are the accepted (expected) options for this marker?
      accepted.opts <- unlist(family_marker.list[[fam]][mname])
      
      # Also allow NA to be an accepted outcome (i.e., it does not indicate an error)
      accepted.opts <- c(accepted.opts, NA) # allow NA in offspring to be accepted
      
      # If the parental genotype was missing, and therefore not able to be evaluated, denote it as such
      if(sum(accepted.opts %in% "missing") > 0){
        
        single_genos_eval.df[,i] <- "parent.missing"
        
        # If the parental genotypes were present, and the offspring genotypes predicted, continue
      }else if(sum(accepted.opts %in% "missing") == 0){
        
        # Are the observed genotypes expected genotypes? 
        single_genos_eval.df[,i] <- single_genos.df[,i] %in% accepted.opts
        
      }
      
      # Then add back offspring NAs as missing
      single_genos_eval.df[is.na(single_genos.df[,i]), i] <- "offspring.NA"
      
    }
    
    # # Confirm the edit for the NAs by checking before and after their removal/addition
    # table(single_genos.df, useNA = "ifany")
    # table(single_genos_eval.df, useNA = "ifany")
    
    write.csv(x = single_genos_eval.df, file = paste0("03_results/offspring_true_and_false_matches_", fam, ".csv"))
    
    
  }
  
  
  #### 04. Prepare final output report ####
  files.FN <- list.files(path = "03_results/", pattern = "offspring_true_and_false_matches_")
  files.FN <- files.FN[grep(pattern = "offspring_true_and_false_matches_all.csv", x = files.FN, invert = T)] # ignore the final output
  files.FN
  
  # Read in all of the T/F match matrices and combine them into one big dataframe
  temp <- NULL; all_data.df <- NULL
  for(i in 1:length(files.FN)){
    
    temp <- read.csv(file = paste0("03_results/", files.FN[i]))
    all_data.df <- rbind(all_data.df, temp)
    
  }
  
  dim(all_data.df)
  all_data.df[1:5,1:10]
  
  # Clean up names a bit
  colnames(all_data.df)[which(colnames(all_data.df)=="X")] <- "indiv"
  colnames(all_data.df) <- gsub(pattern = "^X", replacement = "", x = colnames(all_data.df))
  all_data.df[1:5,1:40]
  #str(all_data.df)
  
  # Improve visualization
  # Convert all to character
  library("dplyr")
  all_data.df <- all_data.df %>%
    mutate_all(as.character)
  all_data.df[1:5,1:10]
  
  # Replace TRUE with dashes to make visually easier to view
  library("tidyverse")
  all_data.df <- all_data.df %>%
    mutate_all(funs(str_replace(., "TRUE", "-")))
  all_data.df[1:5,1:10]
  
  # Remove the first allele indicator
  colnames(all_data.df) <- gsub(pattern = "\\.01", replacement = "", x = colnames(all_data.df))
  all_data.df[1:5,1:10]
  
  
  ## Summarize the results in a dataframe
  # Create matrix with the number of rows being the number of columns (i.e., number of loci evaluated), and four cols
  # Note: the first column is allowed as a dummy
  result.df <- matrix(data = NA, nrow = ncol(all_data.df), ncol = 4)
  dim(result.df) 
  
  # Set column names
  colnames(result.df) <- c("mname", "exp.offsp.geno", "unexp.offsp.geno", "percent.exp.offsp.geno")
  result.df <- as.data.frame(result.df) # make df
  
  # Fill dataframe
  exp.offsp.geno <- NULL; unexp.offsp.geno <- NULL; percent.exp.offsp.geno <- NULL; mname <- NULL
  for(i in 1:ncol(all_data.df)){
    
    mname                 <- colnames(all_data.df)[i]
    result.df[i, "mname"] <- mname
    
    exp.offsp.geno        <- sum(all_data.df[,i]=="-", na.rm = T)
    result.df[i, "exp.offsp.geno"] <- exp.offsp.geno
    
    unexp.offsp.geno      <- sum(all_data.df[,i]=="FALSE", na.rm = T)
    result.df[i, "unexp.offsp.geno"] <- unexp.offsp.geno
    
    percent.exp.offsp.geno      <- exp.offsp.geno / (exp.offsp.geno + unexp.offsp.geno)
    result.df[i, "percent.exp.offsp.geno"] <- percent.exp.offsp.geno
    
  }
  
  # Drop the unneeded column
  result.df <- result.df[result.df$mname!="indiv", ]
  
  # Adjust decimal
  result.df$percent.exp.offsp.geno <- formatC(x = result.df$percent.exp.offsp.geno, digits = 3, format = "f")
  
  head(result.df)
  
  # What is the distribution of the erroneous calls?
  pdf(file = "03_results/hist_per_locus_num_unexpected_genos.pdf", width = 6.5, height = 4)
  hist(result.df$unexp.offsp.geno
       , main = ""
       , las = 1
       , xlab = "Per locus, number of indiv. with unexpected genos"
  )
  dev.off()
  
  
  # How many loci have at least two erroneous calls? 
  # table(result.df$unexp.offsp.geno >= 2) # 95 (pilot), 73 (OCP)
  table(result.df$unexp.offsp.geno >= 4) # 36 (pilot), 35 (OCP)
  
  # Write out per locus info
  write.csv(x = all_data.df, file = paste0("03_results/per_locus_all_results.csv"), row.names = F)
  write.csv(x = result.df, file = paste0("03_results/per_locus_expected_offsp_genos.csv"), row.names = F)
  
  
  # What loci are erroneous in at least two offspring? 
  # loci_to_drop <- result.df[result.df$unexp.offsp.geno >= 2, "mname"]
  loci_to_drop <- result.df[result.df$unexp.offsp.geno >= 4, "mname"]
  
  # write it out
  write.table(x = loci_to_drop, file = "03_results/markers_with_unexpected_genos_in_offspring.txt"
              , quote = F, sep = "\t", row.names = F, col.names = F
  )
  
  
  #### Remove the problematic loci from the original genind file #####
  obj <- obj.bck
  obj
  
  drop_loci(df = obj, drop_file = "03_results/markers_with_unexpected_genos_in_offspring.txt")
  obj_filt
  
  #### Write the new data to a rubias file ####
  # Write out to rubias
  pop_map.FN  # renamed samples
  
  genepop_to_rubias_SNP(data = obj_filt, sample_type = "reference", custom_format = TRUE, micro_stock_code.FN = micro_stock_code.FN
                        , pop_map.FN = pop_map.FN)
  
  print("Your output is available as '03_results/rubias_output_SNP.txt")
  
  
  #### Rename rubias file ####
  # Obtain some variables to create filename
  date <- format(Sys.time(), "%Y-%m-%d")
  indiv.n <- nInd(obj_filt)
  loci.n  <- nLoc(obj_filt)
  rubias_custom.FN <- paste0("03_results/rubias_", indiv.n, "_ind_", loci.n, "_loc_", date, ".txt")
  
  # Copy to rename rubias output file
  file.copy(from = "03_results/rubias_output_SNP.txt", to = rubias_custom.FN, overwrite = T)
  print(paste0("The output is saved as ", rubias_custom.FN))
  
  # Save image
  save.image("03_results/post-filters_prepared_for_parentage_rubias_built_problem_genos_rem_rubias_built.RData")
  
  print("Here you need to copy the above rubias file to amplitools results folder.")
  
  # Go to OCP23_analysis_part_5_2024-02-26.R or 
   

