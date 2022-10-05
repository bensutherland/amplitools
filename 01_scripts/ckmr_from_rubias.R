# Use a rubias file and CKMR to assess parentage
# Sutherland Bioinformatics, 2022-09-19

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Install packages
#install.packages("devtools")
#devtools::install_github("hadley/devtools")
#install.packages("Rcpp")
#install.packages("tidyverse")

require("Rcpp")

### 00. Installing CKMRsim and additional Front Matter
library(devtools)
#devtools::install_github("eriqande/CKMRsim", build_vignettes = TRUE)
library("tidyverse")
library("CKMRsim")

# Load vignettes
vignette("CKMRsim-example-1")

vignette("CKMRsim-example-2-microsatellites")

vignette("CKMRsim-simulating-linked-markers")

vignette("CKMRsim-writing-geno-error-funcs")

# In order to pursue simulations in the face of physical linkage, 
# you must download and install the external dependency, Mendel version 
# 16. For Windows and Mac OS X, CKMRsim will look for the Mendel binary 
# in its default install location. So download it 
# from http://software.genetics.ucla.edu/mendel and do a default 
# install. Note that you need to register with your email in order 
# to download Mendel.


#### Set user variables
input.FN <- "../simple_pop_stats/03_results/rubias_output_SNP.txt"

#### 01. Read in genotype dataset ####
data.df <- read.delim(file = input.FN, header = T, sep = "\t")
data.df[1:10, 1:10]
dim(data.df)

# Only keep parentage samples
data.df <- data.df[data.df$repunit=="VIU_offspring" | data.df$repunit=="VIU_parent",]


# # TODO: Issue: column names apparently cannot be numbers? Not clear where this happened
# test <- read.table(file = input.FN, header = T, sep = "\t")
# test[1:10, 1:10]

# #### SEP COLLS FOR PARENT/ OFFSPRING ####
# # Add column indicating if the sample is offspring or parent
# data.df$hierarchy <- NA
# data.df$hierarchy[grep(pattern = "parent", x = data.df$collection)] <- "parent"
# data.df$hierarchy[grep(pattern = "parent", x = data.df$collection, invert = T)] <- "offspring"
# 
# dim(data.df)
# colnames(data.df)
# data.df[1:100, c("indiv", "hierarchy") ]


# Remove extra cols
data.df <- data.df[, grep(pattern = "sample_type|collection|repunit", x = colnames(data.df), invert = T) ]
dim(data.df)
data.df[1:5, 1:10]


#### SEP COLLS END/ ####
# # Separate parent or offspring datasets
# parent_data.df <- data.df[data.df$hierarchy=="parent",]
# parent_data.df[1:5,1:5]
# dim(parent_data.df)
# parent_data.df <- parent_data.df[, grep(pattern = "hierarchy", x = colnames(parent_data.df), invert = T)]
# dim(parent_data.df)
# 
# offspring_data.df <- data.df[data.df$hierarchy!="parent",]
# offspring_data.df[1:5,1:5]
# dim(offspring_data.df)
# offspring_data.df <- offspring_data.df[, grep(pattern = "hierarchy", x = colnames(offspring_data.df), invert = T)]
# dim(offspring_data.df)
# 
# unique(offspring_data.df[,2])


# For now, keep all data together
genos <- data.df
genos[1:5,1:5]
dim(genos)

#### 02. Compute Allele Frequencies from Genotype Data
# From tutorial
nc <- ncol(genos)
loci <- str_replace(string = names(genos)[seq(2, nc, by = 2)]
                    , pattern = "\\.\\.\\.[0-9]+$"
                    , replacement = ""
                    ) 

length(loci)

# Reset the locus names
names(genos)[seq(2, nc, by = 2)] <- str_c(loci, "1", sep = ".")
names(genos)[seq(3, nc, by = 2)] <- str_c(loci, "2", sep = ".")

genos[1:5,1:5]

# Make into a tibble
genos <- tibble(genos)
genos

# then make some long format genotypes
long_genos <- genos %>% 
  
  gather(key = "loc", value = "Allele", -indiv) %>%
  
  separate(loc, into = c("Locus", "gene_copy"), sep = "\\.") %>%
  
  mutate(Allele = as.character(Allele)) %>%
  
  mutate(Allele = ifelse(Allele == "0", NA, Allele)) %>%
  
  rename(Indiv = indiv)

long_genos

# Generate allele frequencies
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
afreqs_ready <- reindex_markers(alle_freqs)

# Create a CKMR object
kappas

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
PO_U_logls <- extract_logls(ex1_Qs,
                            numer = c(PO = 1),
                            denom = c(U = 1))

PO_U_logls

ggplot(PO_U_logls,
       aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25)


ggplot(PO_U_logls %>% filter(true_relat %in% c("PO", "U")),
       aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25)

FS_U_logls <- extract_logls(ex1_Qs,
                            numer = c(FS = 1),
                            denom = c(U = 1))

ggplot(FS_U_logls %>% filter(true_relat %in% c("FS", "U")),
       aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25)

# Estimating false negative and false positive rates
ex1_PO_is <- mc_sample_simple(ex1_Qs, 
                              nu = "PO",
                              de = "U")

ex1_PO_is

ex1_PO_is_5 <- mc_sample_simple(ex1_Qs, 
                                nu = "PO",
                                de = "U", 
                                lambda_stars = 5)

ex1_PO_is_5

#### note: skipped a few steps #####

# Screen out duplicates
matchers <- find_close_matching_genotypes(LG = long_genos,
                                          CK = ex1_ckmr,
                                          max_mismatch = 6)
matchers

#### Compute Logl Ratios for All Pairwise Comparisons
# Look for parent offspring pairs
indiv_names <- unique(long_genos$Indiv)
parent_ids <- indiv_names[grep(pattern = "BR", x = indiv_names)]
offspring_ids <- indiv_names[grep(pattern = "F", x = indiv_names)]

candidate_parents <- long_genos %>% 
  filter(Indiv %in% parent_ids)

#unique(candidate_parents$Indiv)

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

po_pairwise_logls_greater_than_12 <- po_pairwise_logls %>%
      filter(logl_ratio > 12) %>%
      arrange(desc(logl_ratio))

write.table(x = po_pairwise_logls_greater_than_12, file = "03_results/po_pairwise_logls.txt"
            , sep = "\t", row.names = F, quote = F
            )

fs_pairwise_logls <- pairwise_kin_logl_ratios(D1 = candidate_offspring,
                                              D2 = candidate_offspring,
                                              CK = ex1_ckmr,
                                              numer = "FS",
                                              denom = "U", 
                                              #num_cores = 1
                                              )

fs_pairwise_logls_greater_than_12 <- fs_pairwise_logls %>%
                                    filter(logl_ratio > 12) %>%
                                    arrange(desc(logl_ratio))

write.table(x = fs_pairwise_logls_greater_than_12, file = "03_results/fs_pairwise_logls.txt"
            , sep = "\t", row.names = F, quote = F
)



