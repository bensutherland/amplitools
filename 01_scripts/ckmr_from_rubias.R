# Use a rubias file and CKMR to assess parentage
# Sutherland Bioinformatics, Initialized 2022-09-19

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

## Install packages
#install.packages("devtools")
#devtools::install_github("hadley/devtools")
#install.packages("Rcpp")
#install.packages("tidyverse")

#require("Rcpp")

### 00. Installing CKMRsim and additional Front Matter
library(devtools)
#devtools::install_github("eriqande/CKMRsim", build_vignettes = TRUE)
library("tidyverse")
library("CKMRsim")

# Load vignettes
vignette("CKMRsim-example-1")

#vignette("CKMRsim-example-2-microsatellites")
#vignette("CKMRsim-simulating-linked-markers")
#vignette("CKMRsim-writing-geno-error-funcs")

# ASIDE: re: simulations considering physical linkage #
# In order to pursue simulations in the face of physical linkage, 
# you must download and install the external dependency, Mendel version 
# 16. For Windows and Mac OS X, CKMRsim will look for the Mendel binary 
# in its default install location. So download it 
# from http://software.genetics.ucla.edu/mendel and do a default 
# install. Note that you need to register with your email in order 
# to download Mendel.
# /end/ ASIDE

#### Set user variables
input.FN <- "../simple_pop_stats/03_results/rubias_output_SNP.txt"

#### 01. Read in genotype dataset ####
data.df <- read.delim(file = input.FN, header = T, sep = "\t")
data.df[1:10, 1:10] # note: numeric header receives 'X' prefix
dim(data.df)

# Only keep parentage samples
data.df <- data.df[data.df$repunit=="VIU_offspring" | data.df$repunit=="VIU_parent", ]
data.df[1:10, 1:10]
dim(data.df)

# Remove annotation columns except for the indiv col
data.df <- data.df[, grep(pattern = "sample_type|collection|repunit", x = colnames(data.df), invert = T) ]
data.df[1:5, 1:10]
dim(data.df)

# For now, keep parents and offspring together
genos <- data.df
genos[1:5,1:5]
dim(genos) # indiv, indiv col + nloci x 2
nloci <- (ncol(genos) - 1) / 2

# Reporting
print(paste0("You have ", nrow(genos), " indiv in the dataset"))
print(paste0("You have ", nloci, " loci in the dataset"))


#### 02. Compute Allele Frequencies from Genotype Data ####
nc <- ncol(genos) # note: this will include the indiv ID column first

## Rename columns ##
# Find the names of each locus (only one per allele pair)
loci <- str_replace(string = names(genos)[seq(2, nc, by = 2)]
                    , pattern = "\\.\\.\\.[0-9]+$"
                    , replacement = ""
                    ) 

length(loci)

# Add locus suffix (.1 for allele 1 and .2 for allele 2)
names(genos)[seq(2, nc, by = 2)] <- str_c(loci, "1", sep = ".")
names(genos)[seq(3, nc, by = 2)] <- str_c(loci, "2", sep = ".")

genos[1:5,1:5]

# #### Optional rarefy ####
# # Half loci
# dim(genos)
# genos <- genos[,1:427]
# genos[1:5, 420:427]
# 
# # Quarter loci
# dim(genos)
# genos <- genos[,1:213]
# genos[1:5, 210:213]

# Make into a tibble
genos <- tibble(genos)
genos

# Convert genotypes into long form (from tutorial)
long_genos <- genos %>% 
  
  gather(key = "loc", value = "Allele", -indiv) %>%
  
  separate(loc, into = c("Locus", "gene_copy"), sep = "\\.") %>%
  
  mutate(Allele = as.character(Allele)) %>%
  
  mutate(Allele = ifelse(Allele == "0", NA, Allele)) %>%
  
  rename(Indiv = indiv)

long_genos

# Generate allele frequencies (from tutorial)
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

alle_freqs

# Pass through reindex_markers function to provide a df of AF that CKMRsim will use for sims
# function resets values of AlleIdx and LocIdx in case of removals, rescales AF to sum to 1, and sorts loci along chromosomes
afreqs_ready <- reindex_markers(alle_freqs)

#### 03. Create a CKMR object ####
kappas # data supplied with package

# Assume true and assumed geno error models are the same
# error model is a true-genotype-independent model with error rate (epsilon) equal to 1 in 200
# Distinguish parent-offspring pairs from unrelated, full sib pairs from unrelated, half sibs from full sibs

ex1_ckmr <- create_ckmr(
  D = afreqs_ready,
  kappa_matrix = kappas[c("PO", "FS", "HS", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

ex1_ckmr

# Simulate genotype pairs and calculate log-probabilities
ex1_Qs <- simulate_Qij(ex1_ckmr, 
                       calc_relats = c("PO", "FS", "U"),
                       sim_relats = c("PO", "FS", "HS", "U") )

ex1_Qs

# Compute log-likelihood ratios to view them
# PO, U
PO_U_logls <- extract_logls(ex1_Qs,
                            numer = c(PO = 1),
                            denom = c(U = 1))

PO_U_logls

# Plot densities of logl_ratio for FS, HS, PO, U
pdf(file = "03_results/logl_ratio_u_hs_fs_po.pdf", width = 7, height = 5)
ggplot(PO_U_logls,
       aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25)
dev.off()

# Plot densities for logl_ratio for U, PO (only consider PO and U from true_relat)
pdf(file = "03_results/logl_ratio_u_po.pdf", width = 7, height = 5)
ggplot(PO_U_logls %>% filter(true_relat %in% c("PO", "U")),
       aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25)
dev.off()


# Compute log-likelihood ratios to view them
# FS, U
FS_U_logls <- extract_logls(ex1_Qs,
                            numer = c(FS = 1),
                            denom = c(U = 1))

# Plot densities for logl_ratio for U, FS
pdf(file = "03_results/logl_ratio_u_fs.pdf", width = 7, height = 5)
ggplot(FS_U_logls %>% filter(true_relat %in% c("FS", "U")),
       aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25)
dev.off()

#### 04. Estimate False Positive Rate and False Negative Rate ####
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
# FPR ~3e-14 if we use logl ratio cutoff of 5

# What cutoff do we want to use? 
# How many potential adults and offspring in the study? 
length(grep(pattern = "BR", x = data.df$indiv)) # 63 adults
length(grep(pattern = "BR", x = data.df$indiv, invert = T)) # 102 offspring
63 * 102 # number pairs being tested (6246)
# per pair FPR 3e-14 would leave us with expected number of FP: 
(63 * 102) * 3^-14 # expect ~ 0.001 false positive pairs

# Recommended: require FPR ~10-100 times smaller than reciprocal of the # of comparisons
0.1 * (63 *102) ^ (-1) # 0.00001556 would be a good FPR to aim for, which would mean logl cutoff of 5 should be fine


# View what other logl Lambda_star would provide
ex1_PO_is_5_30 <- mc_sample_simple(ex1_Qs, 
                                   nu = "PO",
                                   de = "U", 
                                   lambda_stars = seq(5, 30, by = 2))
ex1_PO_is_5_30

### It appears, based on the tutorial, that a logl of 5 provides us with a low enough FPR that we could use this
# It also appears that we cannot get the false positive rate any higher than 3.6e-14

#### 05. Estimate FPR and FNR for FS ####
ex1_FS_is <- mc_sample_simple(ex1_Qs, 
                              nu = "FS",
                              de = "U"
                              , lambda_stars = seq(0, 5, by = 0.5)
                              )

ex1_FS_is


#### 06. Empirical detection of PO ####
# Screen out duplicates, for both PO and FS comparisons
matchers <- find_close_matching_genotypes(LG = long_genos,
                                          CK = ex1_ckmr,
                                          max_mismatch = 6)
matchers


#### 07. Compute logl ratios for all pairwise comparisons to look for parent offspring pairs
# Identify parent and offspring IDs
indiv_names <- unique(long_genos$Indiv)
parent_ids <- indiv_names[grep(pattern = "BR", x = indiv_names)]
offspring_ids <- indiv_names[grep(pattern = "BR", x = indiv_names, invert = T)]

candidate_parents <- long_genos %>% 
  filter(Indiv %in% parent_ids)

unique(candidate_parents$Indiv)

candidate_offspring <- long_genos %>% 
  filter(Indiv %in% offspring_ids)

unique(candidate_offspring$Indiv)

po_pairwise_logls <- pairwise_kin_logl_ratios(D1 = candidate_parents, 
                                              D2 = candidate_offspring, 
                                              CK = ex1_ckmr,
                                              numer = "PO",
                                              denom = "U",
                                              #num_cores = 1
                                              )

cutoff <- 5

po_pairwise_logls_over_threshold <- po_pairwise_logls %>%
      filter(logl_ratio > cutoff) %>%
      arrange(desc(logl_ratio))

write.table(x = po_pairwise_logls_over_threshold, file = "03_results/po_pairwise_logls.txt"
            , sep = "\t", row.names = F, quote = F
            )


#### 07. Sibship within offspring #####
fs_pairwise_logls <- pairwise_kin_logl_ratios(D1 = candidate_offspring,
                                              D2 = candidate_offspring,
                                              CK = ex1_ckmr,
                                              numer = "FS",
                                              denom = "U", 
                                              #num_cores = 1
                                              )

fs_pairwise_logls_greater_than_5 <- fs_pairwise_logls %>%
                                    filter(logl_ratio > 5) %>%
                                    arrange(desc(logl_ratio))

write.table(x = fs_pairwise_logls_greater_than_5, file = "03_results/fs_pairwise_logls_greater_than_5.txt"
            , sep = "\t", row.names = F, quote = F
)

# # Observe results closer
# fs.df <- fs_pairwise_logls_greater_than_5
# 
# fs.df <- separate(data = fs.df, col = "D2_indiv"
#          , into = c("fam2", "ind2"), sep = "-", remove = T
#          )
# 
# fs.df <- separate(data = fs.df, col = "D1_indiv"
#                   , into = c("fam1", "ind1"), sep = "-", remove = T
# )
# 
# head(fs.df)
# 
# fs.df$match <- fs.df$fam2==fs.df$fam1
# table(fs.df$match)
# 
# summary(fs.df$logl_ratio[fs.df$match==TRUE])
# summary(fs.df$logl_ratio[fs.df$match==FALSE])
# 
# boxplot(fs.df$logl_ratio ~ fs.df$match
#         , ylab = "logl ratio"
#         , xlab = "Matched family ID")
# 
# # Which families were most likely to be incorrectly matched up?
# incorr.matches <- as.data.frame(fs.df[fs.df$match==FALSE, c("fam2", "fam1")])
# sort(table(paste0(incorr.matches$fam2, "_", incorr.matches$fam1)))


#### 08. Sibship within parents  #####
fs_pairwise_logls_parents <- pairwise_kin_logl_ratios(D1 = candidate_parents,
                                              D2 = candidate_parents,
                                              CK = ex1_ckmr,
                                              numer = "FS",
                                              denom = "U", 
                                              #num_cores = 1
)

fs_pairwise_logls_greater_than_5 <- fs_pairwise_logls_parents %>%
  filter(logl_ratio > 5) %>%
  arrange(desc(logl_ratio))

write.table(x = fs_pairwise_logls_greater_than_5, file = "03_results/fs_pairwise_logls_greater_than_5_parents.txt"
            , sep = "\t", row.names = F, quote = F
)

# # Observe results closer
# fs.df <- fs_pairwise_logls_greater_than_10
# 
# as.data.frame(fs.df)
# 
# summary(fs.df$logl_ratio[fs.df$match==TRUE])
# summary(fs.df$logl_ratio[fs.df$match==FALSE])
# 
# boxplot(fs.df$logl_ratio ~ fs.df$match
#         , ylab = "logl ratio"
#         , xlab = "Matched family ID")
# 
# # Which families were most likely to be incorrectly matched up?
# incorr.matches <- as.data.frame(fs.df[fs.df$match==FALSE, c("fam2", "fam1")])
# sort(table(paste0(incorr.matches$fam2, "_", incorr.matches$fam1)))
# 
# 
# 
