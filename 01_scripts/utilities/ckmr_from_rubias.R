# Assess parentage and estimate power
# Inputs: rubias file
# Sutherland Bioinformatics, Initialized 2022-09-19
# Note: uses a lot of code from the tutorial for CKMRsim, see README for link. 

## requires sourcing 00_initiator.R before starting

# More information: 
#vignette("CKMRsim-example-1")
## ...other vignettes available from CKMRsim page

ckmr_from_rubias <- function(input.FN = "03_prepped_data/cgig_all_rubias.txt", parent_pop = "VIU_F1", offspring_pop = "VIU_F2", parent_pattern = "BR"){
  
  #### 01. Read in genotype dataset ####
  data.df <- read.delim(file = input.FN, header = T, sep = "\t")
  print(data.df[1:10, 1:10]) # note: numeric header receives 'X' prefix
  
  # Reporting
  print(paste0("The data has ", nrow(data.df), " rows and ", ncol(data.df), " columns"))
  
  # Only keep parentage samples
  data.df <- data.df[data.df$repunit==parent_pop | data.df$repunit==offspring_pop, ]
  print(data.df[1:10, 1:10])
  
  # Reporting
  print(paste0("The selected data has ", nrow(data.df), " rows and ", ncol(data.df), " columns"))
  
  # Remove annotation columns except for the indiv col
  print("Removing annotation columns except the one containing individual names")
  data.df <- data.df[, grep(pattern = "sample_type|collection|repunit", x = colnames(data.df), invert = T) ]
  print(data.df[1:5, 1:10])
  
  # Reporting
  print(paste0("The selected, trimmed data has ", nrow(data.df), " rows and ", ncol(data.df), " columns"))
  
  # For now, keep parents and offspring together
  genos <- data.df
  genos[1:5,1:5]
  print(paste0("The data has ", nrow(genos), " individuals and ", (ncol(genos) -1) / 2, " markers"))
  
  
  #### 02. Compute Allele Frequencies from Genotype Data ####
  print("Computing allele frequences")
  nc <- ncol(genos) # note: this will include the indiv ID column first
  
  ## Rename columns ##
  # Find the names of each locus (only one per allele pair)
  loci <- str_replace(string = names(genos)[seq(2, nc, by = 2)]
                      , pattern = "\\.\\.\\.[0-9]+$"
                      , replacement = ""
  ) 
  
  # Add locus suffix (.1 for allele 1 and .2 for allele 2)
  names(genos)[seq(2, nc, by = 2)] <- str_c(loci, "1", sep = ".")
  names(genos)[seq(3, nc, by = 2)] <- str_c(loci, "2", sep = ".")
  
  genos[1:5,1:5]
  
  # ## Optional rarefy ##
  # # Half loci
  # dim(genos)
  # genos <- genos[,1:427]
  # genos[1:5, 420:427]
  # 
  # # Quarter loci
  # dim(genos)
  # genos <- genos[,1:213]
  # genos[1:5, 210:213]
  # ## /END/ Optional rarefy ##
  
  # Make into a tibble
  print("Making genotypes into a tibble")
  genos <- tibble(genos)
  print(genos)
  
  # Convert genotypes into long form (from tutorial)
  print("Converting genotypes to long form")
  long_genos <- genos %>% 
    
    gather(key = "loc", value = "Allele", -indiv) %>%
    
    separate(loc, into = c("Locus", "gene_copy"), sep = "\\.") %>%
    
    mutate(Allele = as.character(Allele)) %>%
    
    mutate(Allele = ifelse(Allele == "0", NA, Allele)) %>%
    
    rename(Indiv = indiv)
  
  print(long_genos)
  
  # Generate allele frequencies (from tutorial)
  print("Generating allele frequencies")
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
  
  print(alle_freqs)
  
  # Pass through reindex_markers function to provide a df of AF that CKMRsim will use for sims
  # function resets values of AlleIdx and LocIdx in case of removals, rescales AF to sum to 1, and sorts loci along chromosomes
  print("Reindexing markers for CKMRsim")
  afreqs_ready <- reindex_markers(alle_freqs)
  
  
  #### 03. Create a CKMR object ####
  print("To create a CKMR object will use the kappas df supplied with CKMRsim")
  print(kappas) # data supplied with package
  
  # Create ckmr object
  print("Creating ckmr object")
  ex1_ckmr <- create_ckmr(
    D = afreqs_ready,
    kappa_matrix = kappas[c("PO", "FS", "HS", "U"), ],
    ge_mod_assumed = ge_model_TGIE,
    ge_mod_true = ge_model_TGIE,
    ge_mod_assumed_pars_list = list(epsilon = 0.005),
    ge_mod_true_pars_list = list(epsilon = 0.005)
  )
  
  print(ex1_ckmr)
  
  # Simulate genotype pairs and calculate log-probabilities
  print("Simulating genotype pairs, then calculating log-probabilities")
  ex1_Qs <- simulate_Qij(ex1_ckmr, 
                         calc_relats = c("PO", "FS", "U"),
                         sim_relats = c("PO", "FS", "HS", "U") )
  
  print(ex1_Qs)
  
  
  #### 04. Compute and extract log-likelihood raios from simulated data
  print("Computing log-likelihood ratios from the simulated data")
  
  # Computing logl: PO, U
  print("Extracting logl from simultated parent-offspring relationship compared to unrelated")
  PO_U_logls <- extract_logls(ex1_Qs,
                              numer = c(PO = 1),
                              denom = c(U = 1))
  
  print(PO_U_logls)
  
  # Plot densities of logl_ratio for FS, HS, PO, U
  print("Plotting densities of logl_ratios for Parent-Offspring and Unrelated")
  print("Saving density plot as 03_results/logl_ratio_u_hs_fs_po.pdf")
  pdf(file = "03_results/logl_ratio_u_hs_fs_po.pdf", width = 7, height = 5)
  ggplot(PO_U_logls,
         aes(x = logl_ratio, fill = true_relat)) +
    geom_density(alpha = 0.25)
  dev.off()
  
  # # Plot densities for logl_ratio for U, PO (only consider PO and U from true_relat)
  # pdf(file = "03_results/logl_ratio_u_po.pdf", width = 7, height = 5)
  # ggplot(PO_U_logls %>% filter(true_relat %in% c("PO", "U")),
  #        aes(x = logl_ratio, fill = true_relat)) +
  #   geom_density(alpha = 0.25)
  # dev.off()
  
  
  # Computing logl: FS, U
  print("Extracting logl from simultated fullsib compared to unrelated")
  FS_U_logls <- extract_logls(ex1_Qs,
                              numer = c(FS = 1),
                              denom = c(U = 1))
  print(FS_U_logls)
  
  # # Plot densities for logl_ratio for U, FS
  # pdf(file = "03_results/logl_ratio_u_fs.pdf", width = 7, height = 5)
  # ggplot(FS_U_logls %>% filter(true_relat %in% c("FS", "U")),
  #        aes(x = logl_ratio, fill = true_relat)) +
  #   geom_density(alpha = 0.25)
  # dev.off()
  
  
  #### 05. Estimate False Positive Rate and False Negative Rate ####
  print("Estimating false positive and false negative rates")
  # PO/U
  
  # Estimating false negative and false positive rates
  # by default will compute FPR assoc. with FNR of 0.3, 0.2, 0.1, 0.05, 0.01, 0.001
  ex1_PO_is <- mc_sample_simple(ex1_Qs, 
                                nu = "PO",
                                de = "U")
  
  ex1_PO_is
  # This shows FPR ~7e-17 when FNR is 0.01
  
  # What would the results look like if we use a logl ratio of 5 as a cutoff (from tutorial)
  ex1_PO_is_5 <- mc_sample_simple(ex1_Qs, 
                                  nu = "PO",
                                  de = "U", 
                                  lambda_stars = 5)
  
  ex1_PO_is_5
  # FPR ~1e-12 if we use logl ratio cutoff of 5
  
  # What cutoff do we want to use? 
  
  ### TODO: needs generalization ###
  
  # How many potential adults and offspring in the study? 
  num_parents <- length(grep(pattern = parent_pattern, x = data.df$indiv))             # 49 adults
  num_offspring <- length(grep(pattern = parent_pattern, x = data.df$indiv, invert = T)) # 102 offspring
  print(paste0("With ", num_parents, " possible parents and ", num_offspring, " possible offspring, "))
  print(paste0("...there are ", num_parents * num_offspring, " pairs being tested.")) # 4998
  
  # per pair FPR 3e-14 would leave us with expected number of FP: 
  print(paste0("Considering the FPR above for a logl ratio of 5, this leaves us with: "))
  print(paste0(    formatC(num_parents * num_offspring * ex1_PO_is_5$FPR[1], format = "e", digits = 2)
               , " possible false positive pairs"))
  
  print("The vignette recommends requiring an FPR that is ~10-100 times smaller than the reciprocal of the number of comparisons, which for this data would be:")
  print(  formatC(0.1 * (num_parents * num_offspring) ^ (-1) , format = "e", digits= 2))
  print("If this value is ~10-100 times larger than your calculated FPR, then you should be good.")
  
  # View what other logl Lambda_star would provide
  print("For comparison, see what other lamda_star values would produce")
  ex1_PO_is_5_30 <- mc_sample_simple(ex1_Qs, 
                                     nu = "PO",
                                     de = "U", 
                                     lambda_stars = seq(5, 30, by = 2))
  print(ex1_PO_is_5_30)
  
  ## Fullsib FPR and FNR
  print("Also estimating FPR and FNR for fullsibs vs. unrelated")
  ex1_FS_is <- mc_sample_simple(ex1_Qs, 
                                nu = "FS",
                                de = "U"
                                , lambda_stars = seq(0, 5, by = 0.5)
  )
  
  print(ex1_FS_is)
  
  # How many potential pairs? 
  ### TODO ###
  
  num_offspring * num_offspring
  
  
  #### 06. Screen out duplicates ####
  matchers <- find_close_matching_genotypes(LG = long_genos,
                                            CK = ex1_ckmr,
                                            max_mismatch = 6)
  print(matchers)
  print("If the above is empty, then there are no close matching individuals (no replicates)")
  
  
  #### 07. Compute logl ratios for all pairwise comparisons to look for parent offspring pairs ####
  # Identify parent and offspring IDs
  indiv_names <- unique(long_genos$Indiv)
  parent_ids <- indiv_names[grep(pattern = parent_pattern, x = indiv_names)]
  offspring_ids <- indiv_names[grep(pattern = parent_pattern, x = indiv_names, invert = T)]
  
  ### NOTE: will need a solution here too, if the indiv name does not hold the brood year
  
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
  print("Estimating pairwise logl ratios from empirical data")
  po_pairwise_logls <- pairwise_kin_logl_ratios(D1 = candidate_parents, 
                                                D2 = candidate_offspring, 
                                                CK = ex1_ckmr,
                                                numer = "PO",
                                                denom = "U",
                                                #num_cores = 1
  )
  
  ### TODO: this should be a flag
  cutoff <- 5
  
  po_pairwise_logls_over_threshold <- po_pairwise_logls %>%
    filter(logl_ratio > cutoff) %>%
    arrange(desc(logl_ratio))
  
  print("Writing out data as '03_results/po_pairwise_logls.txt'")
  write.table(x = po_pairwise_logls_over_threshold, file = "03_results/po_pairwise_logls.txt"
              , sep = "\t", row.names = F, quote = F
  )
  
  
  ## Compute offspring-offspring logls
  fs_pairwise_logls <- pairwise_kin_logl_ratios(D1 = candidate_offspring,
                                                D2 = candidate_offspring,
                                                CK = ex1_ckmr,
                                                numer = "FS",
                                                denom = "U", 
                                                #num_cores = 1
  )
  
  # Keeping only those greater than 5
  fs_pairwise_logls_greater_than_5 <- fs_pairwise_logls %>%
    filter(logl_ratio > cutoff) %>%
    arrange(desc(logl_ratio))
  
  print("Writing out data as '03_results/fs_pairwise_logls_greater_than_5.txt'")
  write.table(x = fs_pairwise_logls_greater_than_5, file = "03_results/fs_pairwise_logls_greater_than_5.txt"
              , sep = "\t", row.names = F, quote = F
  )
  
  
  ## Compute offspring-offspring logls
  fs_pairwise_logls_parents <- pairwise_kin_logl_ratios(D1 = candidate_parents,
                                                        D2 = candidate_parents,
                                                        CK = ex1_ckmr,
                                                        numer = "FS",
                                                        denom = "U", 
                                                        #num_cores = 1
  )
  
  fs_pairwise_logls_parents_greater_than_5 <- fs_pairwise_logls_parents %>%
    filter(logl_ratio > cutoff) %>%
    arrange(desc(logl_ratio))
  
  print("Writing out data as '03_results/fs_pairwise_logls_greater_than_5_parents.txt'")
  write.table(x = fs_pairwise_logls_parents_greater_than_5, file = "03_results/fs_pairwise_logls_greater_than_5_parents.txt"
              , sep = "\t", row.names = F, quote = F
  )
  
  #### TODO: update the names a bit, and ensure FPR and FNR was assessed and reported for PO, FS (offspring), and FS (parents)

}

