# Use the package CKMR-sim to estimate log likelihood distributions and conduct parentage or sibship assignment
#  Largely adapted from CKMR-sim tutorials. Please cite: Anderson EC (2023) CKMRsim. https://github.com/eriqande/CKMRsim.
#  vignette("CKMRsim-example-1")
#  Sutherland Bioinformatics, Initialized 2022-09-19

ckmr_from_rubias <- function(input.FN = "03_prepped_data/cgig_all_rubias.txt", parent_pop = "VIU_F1", offspring_pop = "VIU_F2", cutoff = 5){
  
  # Set variables
  ckmr_results.list <- list()
  
  #### 01. Read in data ####
  # Reporting
  print(paste0("Reading in data from ", input.FN))
  
  # Read in data
  data.df <- read.delim(file = input.FN, header = T, sep = "\t")
  # data.df[1:10, 1:10]
  
  # Reporting
  print(paste0("The data contains ", length(unique(data.df$indiv)), " unique individuals and ", (ncol(data.df)-4)/2, " loci"))
  
  # Retain the number of genotyped loci in each indiv (untyped = NA)
  missing.df <- matrix(data = NA, nrow = nrow(data.df), ncol = 2)
  colnames(missing.df) <- c("indiv", "num_genos")
  missing.df <- as.data.frame(missing.df)
  
  for(row in 1:nrow(data.df)){
    
    missing.df[row,"indiv"] <- data.df[row,"indiv"]
    
    # Note: assumes there are four columns with metadata
    missing.df[row, "num_genos"] <- (ncol(data.df) - sum(is.na(data.df[row,])) - 4)/2
    
  }
  
  # Write out number of genotypes per individual to output
  write.table(x = missing.df, file = "03_results/per_indiv_geno_loci.txt"
              , sep = "\t", row.names = F, col.names = T, quote = F
  )
  
  
  #### 02. Subset dataset to keep only parent and offspring populations
  print(paste0("Keeping only the samples from the parent or offspring populations: ", parent_pop, " and ", offspring_pop))
  data.df <- data.df[data.df$repunit==parent_pop | data.df$repunit==offspring_pop, ]
  
  # Reporting
  print(paste0("The data contains ", length(unique(data.df$indiv)), " unique individuals and ", (ncol(data.df)-4)/2, " loci"))
  
  # Report and export the names of indiv from each generation
  parent_indivs    <- data.df[data.df$repunit==parent_pop, "indiv"]
  offspring_indivs <- data.df[data.df$repunit==offspring_pop, "indiv"]
  print(paste0("The ", length(parent_indivs), " samples in the parental population include: "))
  print(parent_indivs)
  print(paste0("The ", length(offspring_indivs), " samples in the offspring population include: "))
  print(offspring_indivs)
  
  # Export names of analyzed individuals by population
  write.table(x = parent_indivs, file = "03_results/parent_indiv.txt"
              , sep = "\t", row.names = F, col.names = F, quote = F
  )
  
  write.table(x = offspring_indivs, file = "03_results/offspring_indiv.txt"
              , sep = "\t", row.names = F, col.names = F, quote = F
  )
  
  
  #### 03. Restrict to only the indiv col and genotypes ####
  # Remove annotation columns except for the indiv col
  data.df <- data.df[, grep(pattern = "sample_type|collection|repunit", x = colnames(data.df), invert = T) ]
  #print(data.df[1:5, 1:10])
  
  # # Reporting
  # print(paste0("The data contains ", length(unique(data.df$indiv)), " unique individuals and ", (ncol(data.df)-1)/2, " loci"))
  
  
  #### 04. Rename markers to include .1 and .2 suffix ####
  # How many columns in dataset? 
  nc <- ncol(data.df) # note: the col count also includes the 1st col, indiv ID
  
  
  ## Rename marker names in cols ##
  # Find the names of each locus (one per allele pair)
  loci <- str_replace(string = names(data.df)[seq(2, nc, by = 2)]
                      , pattern = "\\.\\.\\.[0-9]+$"
                      , replacement = ""
  ) 
  
  # Add locus suffix (.1 for allele 1 and .2 for allele 2)
  names(data.df)[seq(2, nc, by = 2)] <- str_c(loci, "1", sep = ".")
  names(data.df)[seq(3, nc, by = 2)] <- str_c(loci, "2", sep = ".")
  #print(data.df[1:5,1:5])
  
  
  #### 06. Compute allele frequencies (from CKMRsim tutorial) ####
  print("Computing allele frequences on selected samples")
  
  # Make into a tibble
  data.df <- tibble(data.df)
  #print(data.df)
  
  # Convert genotypes into long form (from CKMRsim tutorial)
  print("Converting genotypes to long form")
  
  long_genos <- data.df %>% 
    
    gather(key = "loc", value = "Allele", -indiv) %>%
    
    separate(loc, into = c("Locus", "gene_copy"), sep = "\\.") %>%
    
    mutate(Allele = as.character(Allele)) %>%
    
    mutate(Allele = ifelse(Allele == "0", NA, Allele)) %>%
    
    rename(Indiv = indiv)
  
  #print(long_genos)
  
  # Calculate allele frequencies (from CKMRsim tutorial)
  alle_freqs <- long_genos %>%
    
    count(Locus, Allele) %>%
    
    group_by(Locus) %>%
    
    mutate(Freq = n / sum(n),
           Chrom = "Unk",
           Pos = as.integer(factor(Locus, levels = loci))) %>%
    
    ungroup() %>%
    
    select(Chrom, Pos, Locus, Allele, Freq) %>%
    
    arrange(Pos, desc(Freq)) %>%
    
    mutate(AlleIdx = NA,
           LocIdx = NA) %>%
    
    filter(!is.na(Allele))
  
  #print(alle_freqs)
  
  # Pass through reindex_markers function to provide a df of AF that CKMRsim will use for sims
  #  function resets values of AlleIdx and LocIdx in case of removals, rescales AF to sum to 1, and sorts loci along chromosomes
  afreqs_ready <- reindex_markers(alle_freqs)
  print("*** Calculated allele frequencies: ")
  print(afreqs_ready)
  
  
  #### 07. Create a CKMR object ####
  print("Create a CKMR object by using the kappas df supplied with CKMRsim")
  #print(kappas) # data supplied with package
  
  # Create ckmr object
  print("Creating ckmr object using kappas PO, FS, HS, and U")
  ex1_ckmr <- create_ckmr(
    D = afreqs_ready,
    kappa_matrix = kappas[c("PO", "FS", "HS", "U"), ],
    ge_mod_assumed = ge_model_TGIE,
    ge_mod_true = ge_model_TGIE,
    ge_mod_assumed_pars_list = list(epsilon = 0.005),
    ge_mod_true_pars_list = list(epsilon = 0.005)
  )
  
  #print(ex1_ckmr)
  
  
  #### 08. Simulate genotype pairs and calculate log-probabilities ####
  print("Simulating genotype pairs")
  ex1_Qs <- simulate_Qij(ex1_ckmr, 
                         calc_relats = c("PO", "FS", "U"),
                         sim_relats = c("PO", "FS", "HS", "U") )
  
  # Compute and extract log-likelihood raios from simulated data
  print("Computing log-likelihood ratios from the simulated data")
  
  
  ##### 08.1 PO vs. U #####
  # Computing logl: PO, U
  print("Simulated parent-offspring (PO) vs. Unrelated (U)")
  PO_U_logls <- extract_logls(ex1_Qs,
                              numer = c(PO = 1),
                              denom = c(U = 1))
  
  #print(PO_U_logls)
  
  # Plot densities of logl_ratio for PO, U
  po_logl_density_plot.FN <- "03_results/logl_ratio_PO_U.pdf"
  print(paste0("Plotting densities of logl_ratios, saving to ", po_logl_density_plot.FN))
  
  p <- ggplot(PO_U_logls,
         aes(x = logl_ratio, fill = true_relat)) +
       geom_density(alpha = 0.25)
  
  pdf(file = po_logl_density_plot.FN, width = 7, height = 5)
  print(p)
  dev.off()
  
  # # Plot densities for logl_ratio for U, PO (only consider PO and U from true_relat)
  # pdf(file = "03_results/logl_ratio_u_po.pdf", width = 7, height = 5)
  # ggplot(PO_U_logls %>% filter(true_relat %in% c("PO", "U")),
  #        aes(x = logl_ratio, fill = true_relat)) +
  #   geom_density(alpha = 0.25)
  # dev.off()
  
  
  ##### 08.2 FS vs. U #####
  # Computing logl: FS, U
  print("Simulated full-sib (PO) vs. unrelated (U)")
  FS_U_logls <- extract_logls(ex1_Qs,
                              numer = c(FS = 1),
                              denom = c(U = 1))
  #print(FS_U_logls)
  
  
  # Plot densities of logl_ratio for FS, U
  fs_logl_density_plot.FN <- "03_results/logl_ratio_FS_U.pdf"
  print(paste0("Plotting densities of logl_ratios, saving to ", fs_logl_density_plot.FN))
  
  p <- ggplot(FS_U_logls,
              aes(x = logl_ratio, fill = true_relat)) +
    geom_density(alpha = 0.25)
  
  pdf(file = fs_logl_density_plot.FN, width = 7, height = 5)
  print(p)
  dev.off()
  
  # # Plot densities for logl_ratio for U, FS
  # pdf(file = "03_results/logl_ratio_u_fs.pdf", width = 7, height = 5)
  # ggplot(FS_U_logls %>% filter(true_relat %in% c("FS", "U")),
  #        aes(x = logl_ratio, fill = true_relat)) +
  #   geom_density(alpha = 0.25)
  # dev.off()
  
  
  #### 09. Estimate false positive rates and false negative rates ####
  # Retain info
  ckmr_results.list[["cutoff_applied"]]    <- cutoff
  
  ##### 09.1 Parent-offspring comparisons #####
  print("Estimating false positive and false negative rates for parent-offspring comparisons")
  
  # Evaluate expected FPR/FNR with default FNR range and user-set cutoff (default FNR: 0.01, 0.05, 0.1, 0.2, 0.3)
  ex1_PO_is_5 <- mc_sample_simple(ex1_Qs, 
                                  nu = "PO",
                                  de = "U", 
                                  lambda_stars = cutoff)
  print(ex1_PO_is_5)
  
  # Retain info
  ckmr_results.list[["PO_FPRs"]]           <- ex1_PO_is_5
  ckmr_results.list[["PO_FPR_set_cutoff"]] <- formatC(x = ex1_PO_is_5$FPR[1], format = "e", digits = 2)
  print(paste0("The FPR at the set cutoff is ** ", ckmr_results.list[["PO_FPR_set_cutoff"]], " **"))
  
  # Summarize and record the number of comparisons applied
  num_parents   <- length(parent_indivs)
  num_offspring <- length(offspring_indivs)
  ckmr_results.list[["PO_number_pairs_tested"]] <- num_parents * num_offspring
  print(paste0("With ", num_offspring, " offspring and ", num_parents
               , " parents, there are ", ckmr_results.list[["PO_number_pairs_tested"]]
               , " pairs being tested."))
  
  # How many false positives are expected in the dataset with the chosen cutoff? 
  FP_expected <- ckmr_results.list[["PO_number_pairs_tested"]] * as.numeric(ckmr_results.list[["PO_FPR_set_cutoff"]])
  ckmr_results.list[["PO_estim_FP_pairs"]] <- formatC(x = FP_expected, format = "e", digits = 2)
  print(paste0("The expected number of false positives in the dataset is: ", ckmr_results.list[["PO_estim_FP_pairs"]]))
  
  # The vignette recommends considering the reciprocal of the number of comparisons being made
  recip_of_comps <- 0.1 * ckmr_results.list[["PO_number_pairs_tested"]] ^ (-1)
  ckmr_results.list[["PO_reciproc_num_comps"]] <- formatC(x = recip_of_comps, format = "e", digits = 2)
  
  # Compare between the two
  print("The calculated FPR at user-set logl cutoff: ")
  print(ckmr_results.list[["PO_FPR_set_cutoff"]])
  print("0.1x the reciprocal of the number of comparisons is: ")
  print(ckmr_results.list[["PO_reciproc_num_comps"]])
  print("If the calc FPR is smaller than the reciprocal number of comparisons, then proceed.")
  
  # Also retain other lambda star values
  print("For comparison, also see what other lambda_star values would produce: ")
  ex1_PO_is_5_30 <- mc_sample_simple(ex1_Qs, 
                                     nu = "PO",
                                     de = "U", 
                                     lambda_stars = seq(cutoff, 30, by = 2))
  
  ckmr_results.list[["PO_range_FNR_FPR"]] <- ex1_PO_is_5_30
  print(ckmr_results.list[["PO_range_FNR_FPR"]])
  
  
  ##### 09.2 Full-sib comparisons #####
  # Reporting
  print("Estimating false positive and false negative rates for full-sib comparisons")
  
  ex1_FS_is <- mc_sample_simple(ex1_Qs, 
                                nu = "FS",
                                de = "U"
                                , lambda_stars = cutoff
  )
  
  print(ex1_FS_is)
  
  # Retain info
  ckmr_results.list[["FS_FPRs"]]           <- ex1_FS_is
  ckmr_results.list[["FS_FPR_set_cutoff"]] <- formatC(x = ex1_FS_is$FPR[1], format = "e", digits = 2)
  print(paste0("The FPR at the set cutoff is ** ", ckmr_results.list[["FS_FPR_set_cutoff"]], " **"))
  
  # How many potential pairs? 
  print(paste0("With ", num_offspring, " offspring, there are ", num_offspring * num_offspring, " pairs being tested."))
  print("100x smaller than the reciprocal of the number of pairs = ")
  print(  formatC(0.1 * (num_offspring * num_offspring) ^ (-1) , format = "e", digits= 2))
  print("If your FPR is smaller than the above, then proceed.")
  
  print("***Completed selection of logl cutoffs***")
  

  #### 10. Screen out duplicates ####
  print("Check for dupl. indiv. that exist in the dataset using the 'find_close_matching_genotypes() function'")
  matchers <- find_close_matching_genotypes(LG = long_genos,
                                            CK = ex1_ckmr,
                                            max_mismatch = 6
                                            )
  print(matchers)
  print("If the above printed df is empty, then there are no close matching individuals (no replicate individuals)")
  ## TODO: how to screen out if they are present? ##
  
  
  #### 11. Compute logl ratios for all pairwise comparisons to look for parent offspring pairs ####
  print("Computing logl ratios for all pairwise comparisons to identify parent-offspring pairs")

  # Identify parent and offspring IDs
  indiv_names   <- unique(long_genos$Indiv)
  parent_ids    <- indiv_names[indiv_names %in% parent_indivs] # do this way in case filtering happened between identifying IDs and now (e.g., duplicates)
  offspring_ids <- indiv_names[indiv_names %in% offspring_indivs]
  
  # Clean space
  rm(parent_indivs)
  rm(offspring_indivs)
  
  print("Separating parent and offspring data")
  
  # Isolate parent data
  candidate_parents <- long_genos %>% 
    filter(Indiv %in% parent_ids)
  
  print(paste0("These are your unique parents: "))
  print(unique(candidate_parents$Indiv))
  
  # Isolate offspring data
  candidate_offspring <- long_genos %>% 
    filter(Indiv %in% offspring_ids)
  
  print(paste0("These are your unique offspring: "))
  print(unique(candidate_offspring$Indiv))
  
  # Compute parent-offspring logls
  print("Estimating PO pairwise logl ratios from empirical data")
  po_pairwise_logls <- pairwise_kin_logl_ratios(D1 = candidate_parents, 
                                                D2 = candidate_offspring, 
                                                CK = ex1_ckmr,
                                                numer = "PO",
                                                denom = "U",
                                                #num_cores = 1
  )
  

  # Keep pairwise that are above the logl threshold
  print(paste0("Keeping pairwise logls above set threshold logl = ", cutoff))
  po_pairwise_logls_over_threshold <- po_pairwise_logls %>%
    filter(logl_ratio > cutoff) %>%
    arrange(desc(logl_ratio))
  
  # Write out results
  po_output.FN <- paste0("03_results/po_", parent_pop, "_vs_", offspring_pop, "_pw_logl_", cutoff, ".txt")
  print(paste0("Writing out data as ", po_output.FN))
  write.table(x = po_pairwise_logls_over_threshold, file = po_output.FN
              , sep = "\t", row.names = F, quote = F
  )
  
  # Write out all raw results (including below threshold)
  print("Also writing out all results")
  write.table(x = po_pairwise_logls, file = "03_results/po_pairwise_all_no_cutoff.txt"
              , sep = "\t", row.names = F, quote = F
              )
  
  ## Compute offspring-offspring logls
  print("Estimating FS (offspring) pairwise logl ratios from empirical data")
  fs_pairwise_logls <- pairwise_kin_logl_ratios(D1 = candidate_offspring,
                                                D2 = candidate_offspring,
                                                CK = ex1_ckmr,
                                                numer = "FS",
                                                denom = "U", 
                                                #num_cores = 1
  )
  
  # Keeping only those greater than cutoff
  fs_pairwise_logls_over_threshold <- fs_pairwise_logls %>%
    filter(logl_ratio > cutoff) %>%
    arrange(desc(logl_ratio))
  
  # Write out results
  fs_output.FN <- paste0("03_results/fs_offsp_", offspring_pop, "_pw_logl_", cutoff, ".txt")
  
  print(paste0("Writing out data as ", fs_output.FN))
  write.table(x = fs_pairwise_logls_over_threshold, file = fs_output.FN
              , sep = "\t", row.names = F, quote = F
  )
  
  
  ## Compute parental FS logls
  print("Estimating FS (parent) pairwise logl ratios from empirical data")
  fs_pairwise_logls_parents <- pairwise_kin_logl_ratios(D1 = candidate_parents,
                                                        D2 = candidate_parents,
                                                        CK = ex1_ckmr,
                                                        numer = "FS",
                                                        denom = "U", 
                                                        #num_cores = 1
  )
  
  # Keeping only those greater than cutoff
  fs_pairwise_logls_parents_over_threshold <- fs_pairwise_logls_parents %>%
    filter(logl_ratio > cutoff) %>%
    arrange(desc(logl_ratio))
  
  # Write out results
  fs_output_parents.FN <- paste0("03_results/fs_parent_", parent_pop, "_pw_logl_", cutoff, ".txt")
  print(paste0("Writing out data as ", fs_output_parents.FN))
  write.table(x = fs_pairwise_logls_parents_over_threshold, file = fs_output_parents.FN
              , sep = "\t", row.names = F, quote = F
  )
  
  # Write out ckmr_results.list
  output.FN <- paste0("03_results/ckmr_run_logl_", cutoff, "_result_log.txt")
  capture.output(ckmr_results.list, file = output.FN)
  
  ##### 12. Prepare final report ####
  prep_report(relationship = "PO", input.FN = po_output.FN, offspring_ids = offspring_ids)

}

