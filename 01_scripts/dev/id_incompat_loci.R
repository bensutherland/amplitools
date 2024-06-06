# Identify Mendelian incompatibility loci based on empirical data
# 2024-06-07; Ben J. G. Sutherland (VIU, SBIO)
# Requires that 01_scripts/utilities/ckmr_from_rubias.R has been run
# Assumes diploid, assumes no phasing

id_incompat_loci <- function(data = obj, input_report.FN = "03_results/ckmr_input_rubias_mhaps_396_2024-05-24_VIU_F2_vs_VIU_F1_2024-06-05/po_VIU_F1_vs_VIU_F2_pw_logl_5_report.txt"
                             , strong_logl_cutoff = 10, strong_parent_obs = 2, num_offsp_incompat = 4){
 
  #### Part 1. Identify strong empirical trios as 'true' relationships ####
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
  
  
  #### Part 2. Determine expected offspring genos ####
  #### TODO: this should be a separate function ####
  # What are the unique parent combinations? 
  head(empirical_family_map.df)
  dim(empirical_family_map.df)
  
  # Create a df containing the unique strong trio parental pairs
  unique_families.df <- empirical_family_map.df[!duplicated(x = empirical_family_map.df$family.id), ]
  unique_families.df <- unique_families.df[, c("family.id", "parent1", "parent2")]
  head(unique_families.df)
  
  # Create empty matrix to fill with expected possible offspring genos based on parents
  loci <- matrix(data = NA, nrow = nrow(unique_families.df), ncol = nLoc(data_strong))
  colnames(loci) <- locNames(data_strong)
  loci[1:5,1:5]
  unique_fams_and_loci.df <- cbind(unique_families.df, loci)
  unique_fams_and_loci.df <- as.data.frame(unique_fams_and_loci.df)
  unique_fams_and_loci.df[1:5,1:5] # this will be filled
  dim(unique_fams_and_loci.df)
  
  ## Determine the potential offspring genotypes per pairing
  # set nulls
  p1 <- NULL; p2 <- NULL; fam <- NULL; family_marker.list <- list()
  for(i in 1:nrow(unique_fams_and_loci.df)){
    
    # Set variables
    fam <- as.character(unique_fams_and_loci.df$family.id[i])
    p1 <-  as.character(unique_fams_and_loci.df$parent1[i])
    p2 <-  as.character(unique_fams_and_loci.df$parent2[i])
    
    # Reporting
    print(paste0("Analyzing empirical family ", fam, ", comprised of parents: ", p1, " and ", p2))
    
    loi <- NULL; p1.alleles <- NULL; p2.alleles <- NULL; offspring_genos <- NULL
    opt1 <- NULL; opt2 <- NULL; opt3 <- NULL; opt4 <- NULL
    for(l in 1:nLoc(data_strong)){
      
      # Set locus of interest
      loi <- locNames(data_strong)[l]
    
      # Note: how to determine if either parent has missing data for the locus
      #!isTRUE(all(is.na(data_strong[p1, loc = loi]$tab))) # present
      #!isTRUE(all(is.na(data_strong[p2, loc = loi]$tab))) # missing
      
      # If both parents have complete data for the locus, determine potential offspring genotypes
      #  note: assumes that if an allele is NA, all alleles for that locus are NA (hence the '1' below)
      if(!is.na(data_strong[p1, loc = loi]$tab[1]) & !is.na(data_strong[p2, loc = loi]$tab[1])){
        
        # Obtain the parental alleles that have a non-zero status (present as 1 or 2)
        p1.alleles <- colnames(data_strong[p1, loc = loi]$tab)[which(data_strong[p1, loc = loi]$tab > 0)]
        p2.alleles <- colnames(data_strong[p2, loc = loi]$tab)[which(data_strong[p2, loc = loi]$tab > 0)]
        
        # if either parent shows only a single allele for locus, assume locus is homozygous
        if(length(p1.alleles)==1){
          p1.alleles <- c(p1.alleles, p1.alleles)
        }
        
        if(length(p2.alleles)==1){
          p2.alleles <- c(p2.alleles, p2.alleles)
        }
        
        # Determine putative offspring genotypes? 
        #  note: alleles are sorted for matching purposes
        #### TODO: is the sorting OK? ####
        opt1 <- sort(c(p1.alleles[1], p2.alleles[1]))
        opt1 <- paste0(opt1, collapse = ";")
        
        opt2 <- sort(c(p1.alleles[1], p2.alleles[2]))
        opt2 <- paste0(opt2, collapse = ";")
        
        opt3 <- sort(c(p1.alleles[2], p2.alleles[1]))
        opt3 <- paste0(opt3, collapse = ";")
        
        opt4 <- sort(c(p1.alleles[2], p2.alleles[2]))
        opt4 <- paste0(opt4, collapse = ";")
        
        offspring_genos <- paste0(opt1, ";;", opt2, ";;", opt3, ";;", opt4)
        
        # Save offspring genos into the matrix at the appropriate location
        unique_fams_and_loci.df[unique_fams_and_loci.df$family.id==fam, loi] <- offspring_genos
        
       
        # if either parent is missing, set the potential offspring genotypes as missing 
      }else{
        
        # set this as missing data
        offspring_genos <- paste0(NA, ";;", NA, ";;", NA, ";;", NA)
        
        unique_fams_and_loci.df[unique_fams_and_loci.df$family.id==fam, loi] <- offspring_genos
        
        
      }
      
      
    }
  }
    
  # Check output
  unique_fams_and_loci.df[1:5, 1:5]
   
  
  #### THIS COULD BE A STOP POINT FOR A FUNCTION ####
  
  # Save out parental genos per marker
  write.csv(x = unique_fams_and_loci.df, file = "03_results/parental_genos_per_marker.csv"
            , quote = F, row.names = F
            )
  #### TODO: rename this, its not the parents, its the actual compatible genos ####
    
    

  #### Part 3. Evaluate offspring against parents ####
  #### TODO: this should be a separate function ####
    
  # Who are the offspring of the family? 
  #indNames(data_strong) # individuals available
  head(empirical_family_map.df)
  
  # What are the unique families in the strong trios? 
  families <- unique(empirical_family_map.df$family.id)
  
  # Build an empty matrix for scored offspring genotypes (correct or not based on parents)
  offspring_genos.mat <- matrix(data = NA, nrow = nrow(empirical_family_map.df), ncol = nLoc(data_strong))
  colnames(offspring_genos.mat) <- locNames(data_strong)
  empirical_family_map_w_scores.df <- cbind(empirical_family_map.df, offspring_genos.mat)
  empirical_family_map_w_scores.df[1:5,1:5]
  
  # Build an empty matrix to also save the obs offspring geno
  empirical_family_map_w_offspring_genos.df <- empirical_family_map_w_scores.df
  empirical_family_map_w_offspring_genos.df[1:5,1:5]
  
  
  fam <- NULL; p1 <- NULL; p2 <- NULL; keep <- NULL; offspring <- NULL
  for(o in 1:length(families)){
    
    # Set family
    fam <- families[o]
    p1 <- as.character(gsub(pattern = "__.*", replacement = "", x = fam))
    p2 <- as.character(gsub(pattern = ".*__", replacement = "", x = fam))
    
    # Reporting
    print(paste0("Analyzing family ", fam, ", comprised of parents: ", p1, " and ", p2))
    print(paste0("Offspring to compare are based on the strong trios dataframe"))
    
    # Keep strong trio offspring of this family
    keep <- empirical_family_map.df[empirical_family_map.df$family.id==fam, "offspring"]
    offspring <- data_strong[keep]
    print(indNames(offspring))
    
    # Loop over to check the offspring's results
    ooi <- NULL
    for(i in 1:nInd(offspring)){
      
      # Define the offspring
      ooi <- indNames(offspring)[i]
      print(ooi)
      
      # Define the family
      # foi <- empirical_family_map.df[empirical_family_map.df$offspring==ooi, "family.id"] # technically not necessary, we define above
      
      loi <- NULL; offspring_geno <- NULL; opts <- NULL
      for(l in 1:nLoc(offspring)){
        
        # Define the locus
        loi <- locNames(offspring)[l]
        
        # What is the offspring's geno for this locus?
        offspring_geno <- colnames(offspring[ooi, loc=loi]$tab)[which(offspring[ooi, loc=loi]$tab > 0)]
        
        # If homozygous, it will be length 1, so duplicate the allele to make diploid
        if(length(offspring_geno)==1){
          
          offspring_geno <- c(offspring_geno, offspring_geno)
          
        }
        
        # If offspring geno is missing, set as such and save to the scores file
        if(length(offspring_geno)==0){
        
          offspring_geno <- "offspring_missing"
          
          # and save
          empirical_family_map_w_scores.df[empirical_family_map_w_scores.df$offspring==ooi, loi] <- offspring_geno
         
        # If offspring geno is present, evaluate and save
        }else{
          
          # Sort for matching purposes
          offspring_geno <- paste0(sort(offspring_geno), collapse = ";")
          
          # Determine if it is one of the acceptable allele combinations from the parental dataset
          opts <- strsplit(x = unique_fams_and_loci.df[unique_fams_and_loci.df$family.id==fam, loi], split = ";;")
          opts <- unlist(opts)
          opts
          
          # If the locus is missing in the parents, set as missing here
          if(all((unlist(strsplit(x = opts, split = ";;"))=="NA"))){
            
            empirical_family_map_w_scores.df[empirical_family_map_w_scores.df$offspring==ooi, loi] <- "parent_missing"
            
          # If the locus is present in parents, score
          }else{
            
            # if the allelic combination is present in parental opts, give it a TRUE
            if(offspring_geno %in% opts){
              
              empirical_family_map_w_scores.df[empirical_family_map_w_scores.df$offspring==ooi, loi] <- "OK"
              
              # if the allelic combination is not in parental opts, save the geno
            }else if(!offspring_geno %in% opts){
              
              empirical_family_map_w_scores.df[empirical_family_map_w_scores.df$offspring==ooi, loi] <- offspring_geno
              
            }
          }
          
        }

        # Aside: save the geno in the geno-storing df
        empirical_family_map_w_offspring_genos.df[empirical_family_map_w_offspring_genos.df$offspring==ooi, loi] <- offspring_geno
        
      }
    }
  }

    
  
  empirical_family_map_w_scores.df[1:8,1:10]
  empirical_family_map_w_offspring_genos.df[1:5,1:10]
  

  # Save out
  write.csv(x = empirical_family_map_w_scores.df, file = "03_results/offspring_compats_per_marker_scores.csv"
            , quote = F, row.names = F
  )
  
  # Save out
  write.csv(x = empirical_family_map_w_offspring_genos.df, file = "03_results/offspring_compats_per_marker_genos.csv"
            , quote = F, row.names = F
  )
  
  
  
  
  
  #### Tally outcome ####
  per_locus_tally.df <- as.data.frame(locNames(data_strong))
  # colnames(empirical_family_map_w_scores.df)[grep(pattern = "offspring|family.id|parent1|parent2", x = colnames(empirical_family_map_w_scores.df), invert = T)]
  colnames(per_locus_tally.df) <- "locus"
  head(per_locus_tally.df)
  
  per_locus_tally.df$number_OK <- NA
  per_locus_tally.df$offspring_missing <- NA
  per_locus_tally.df$parent_missing <- NA
  per_locus_tally.df$number_incompats <- NA
  per_locus_tally.df$percent_incompat <- NA
  
  loi <- NULL; number_OK <- NULL; offspring_missing <- NULL; parent_missing <- NULL; total_cases <- NULL; number_incompats <- NULL; percent_incompat <- NULL
  for(i in 1:nrow(per_locus_tally.df)){
    
    loi <- per_locus_tally.df[i,"locus"]
    
    number_OK          <- sum(empirical_family_map_w_scores.df[,loi]=="OK", na.rm = T)
    offspring_missing  <- sum(empirical_family_map_w_scores.df[,loi]=="offspring_missing", na.rm =T)
    parent_missing     <- sum(empirical_family_map_w_scores.df[,loi]=="parent_missing", na.rm = T)
    total_cases        <- length(empirical_family_map_w_scores.df[,loi])
    
    number_incompats <- total_cases - (number_OK + offspring_missing + parent_missing)
    percent_incompat <- number_incompats / (number_OK+number_incompats)
    
    
    # Save
    per_locus_tally.df[per_locus_tally.df$locus==loi, "number_OK"] <- number_OK
    per_locus_tally.df[per_locus_tally.df$locus==loi, "offspring_missing"]  <- offspring_missing
    per_locus_tally.df[per_locus_tally.df$locus==loi, "parent_missing"] <- parent_missing
    per_locus_tally.df[per_locus_tally.df$locus==loi, "total_cases"]    <- total_cases
    per_locus_tally.df[per_locus_tally.df$locus==loi, "number_incompats"]    <- number_incompats
    per_locus_tally.df[per_locus_tally.df$locus==loi, "percent_incompat"]    <- percent_incompat
    
  }
  
  head(per_locus_tally.df) 
  
  
  hist(per_locus_tally.df$number_incompats, breaks = 20)
  sum(per_locus_tally.df$number_incompats > 4, na.rm = T)
  nrow(per_locus_tally.df)
  
  nrow(per_locus_tally.df) -   sum(per_locus_tally.df$number_incompats > 4, na.rm = T)
  # 223 loci remain
  

  write.csv(x = per_locus_tally.df, file = "03_results/per_locus_tally_incompats.csv", quote = F
            , row.names = F
            )
  
  loci_to_drop <- per_locus_tally.df[per_locus_tally.df$number_incompats > 4, "locus"]
  length(loci_to_drop)
  write.table(x = loci_to_drop, file = "03_results/incompat_loci.csv", quote = F, sep = ","
            , col.names = F, row.names = F)
  
  
}
    
  #### FRONT EDGE ####  
    

  
  
  
  
  
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
   

