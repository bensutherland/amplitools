# Script to prepare submission form for the rhAmp assay
# B. Sutherland (VIU, SBIO)

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("tidyr")
library("tidyr")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts\\/designer", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# User set variables
species <- "cgig"

# Global variables (do not change)
resources.dir <- "02_input_data"
input.dir     <- "10_designer"
output.dir    <- "10_designer"

# Does the marker name have full information about position etc (RADseq) or is it just a marker name? 
full_info <- FALSE # just a marker name


#### 1. Import data ####
### Import sequence data ###
seq.df <- read.delim2(file = paste0(input.dir, "/selected_chr_and_seq.txt"), header = F)
seq.df <- as.data.frame(x = seq.df)
colnames(seq.df) <- c("chr_info", "seq")
head(seq.df, n = 1)
tail(seq.df, n = 3)

# Separate columns into individual details
seq.df <- separate(data = seq.df, col = "chr_info", into = c("mname", "chr_info"), sep = "::", remove = T)
head(seq.df)

# Substitute the generic mname with the chosen species name
seq.df$mname <- gsub(pattern = "mname", replacement = species, x = seq.df$mname)
head(seq.df)

# Separate columns into additional individual details
seq.df <- separate(data = seq.df, col = "chr_info", into = c("chr", "pos_range"), sep = ":", remove = T)
seq.df <- separate(data = seq.df, col = "pos_range", into = c("lower_range", "upper_range"), sep = "-", remove = T)

head(x = seq.df, n = 2)
tail(x = seq.df, n = 3)

# Update formats
seq.df$seq         <- as.character(seq.df$seq)
seq.df$lower_range <- as.numeric(seq.df$lower_range)
seq.df$upper_range <- as.numeric(seq.df$upper_range)
head(x = seq.df, n = 1)
str(seq.df)


### Import marker data ###
# This marker info was taken from the VCF, that was taken from the original genotyping
minfo.df <- read.delim2(file = paste0(input.dir, "/vcf_selection.csv"), header = F, sep = ",")
minfo.df <- as.data.frame(x = minfo.df, stringsAsFactors = F)
head(minfo.df)
tail(minfo.df)
colnames(x = minfo.df) <- c("chr", "pos", "info", "ref", "alt")

minfo.df$ref <- as.character(minfo.df$ref)
minfo.df$alt <- as.character(minfo.df$alt)
minfo.df$chr <- as.character(minfo.df$chr)

head(minfo.df, n = 3)
tail(minfo.df, n = 3)
str(minfo.df)

# Separate columns into individual details
if(full_info == TRUE){
  
  minfo.df <- separate(data = minfo.df, col = "info", into = c("mname", "SNP_pos", "align_strand"), sep = ":", remove = T)
  
}else if(full_info == FALSE){
  
  minfo.df$mname <- minfo.df$info
  
}


head(minfo.df, n = 2)
tail(minfo.df, n = 3)
str(minfo.df)

# Add species name to mname
minfo.df$mname <- paste0(species, "_", minfo.df$mname)
head(minfo.df, n = 2)


#### 02. Combine datasets into one ####
## Merge data
head(x = seq.df, n = 1)
tail(x = seq.df, n = 3)
head(minfo.df, n = 1)
nrow(seq.df)
nrow(minfo.df)

seq_and_minfo.df <- merge(x = seq.df, y = minfo.df, by = "mname")
nrow(seq_and_minfo.df)

head(seq_and_minfo.df, n = 2)

# Change nucleotides to all upper case
seq_and_minfo.df$seq <- toupper(x = seq_and_minfo.df$seq)

# Separate the sequence into the flanking 5' 200 bp, the target (genome's ref allele), then the flanking 3' 200 bp
seq_and_minfo.df$left_seq  <- substr(x = seq_and_minfo.df$seq, start =   1, stop = 200) 
seq_and_minfo.df$ref_nuc   <- substr(x = seq_and_minfo.df$seq, start = 201, stop = 201)
seq_and_minfo.df$right_seq <- substr(x = seq_and_minfo.df$seq, start = 202, stop = 401)

# Data checking
head(seq_and_minfo.df, n = 1)
tail(seq_and_minfo.df, n = 3)
seq_and_minfo.df[1:2, c("ref", "alt", "ref_nuc")]

# How many didn't match the expectation for either the ref or alt nucl?
table(seq_and_minfo.df$ref_nuc == seq_and_minfo.df$ref | seq_and_minfo.df$ref_nuc == seq_and_minfo.df$alt)
# That should be all TRUE
seq_and_minfo.df[which((seq_and_minfo.df$ref_nuc == seq_and_minfo.df$ref | seq_and_minfo.df$ref_nuc == seq_and_minfo.df$alt)==FALSE), ]
# There should be no output to this


# Re-combine the seq data back together
seq_and_minfo.df$formatted_seq <- paste0(seq_and_minfo.df$left_seq
                                         , "[", seq_and_minfo.df$ref, "/", seq_and_minfo.df$alt, "]"
                                         , seq_and_minfo.df$right_seq
)


# Write out results (full) for troubleshooting
write.csv(x = seq_and_minfo.df, file = paste0(output.dir, "/seq_and_minfo_all_data.csv"), quote = F, row.names = F)

# rhAmp only needs name and seq with alleles
seq_and_minfo.df <- seq_and_minfo.df[, c("mname", "formatted_seq")]

head(seq_and_minfo.df)

write.csv(x = seq_and_minfo.df, file = paste0(output.dir, "/seq_and_minfo_for_submission.csv"), quote = F, row.names = F)



#### Next Steps: #####
# Can do some data checking, and then can submit the markers for designing primers
