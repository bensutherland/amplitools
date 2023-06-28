# Script to prepare submission form for amplicon panel
# B. Sutherland (Sutherland Bioinformatics)

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
species <- "myes"

# Global variables (do not change)
resources.dir <- "02_input_data"
input.dir     <- "10_designer"
output.dir    <- "10_designer"


#### 1. Import data ####
### Import sequence data ###
seq.df <- read.delim2(file = paste0(input.dir, "/selected_chr_and_seq.txt"), header = F)
seq.df <- as.data.frame(x = seq.df)
colnames(seq.df) <- c("chr_info", "seq")
head(seq.df, n = 1)
tail(seq.df, n = 3)

# Separate columns into individual details
seq.df <- separate(data = seq.df, col = "chr_info", into = c("mname", "chr_info"), sep = "::", remove = T)
seq.df$mname <- gsub(pattern = "mname", replacement = species, x = seq.df$mname)
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

head(minfo.df, n = 2)
tail(minfo.df, n = 3)
str(minfo.df)

# Separate columns into individual details
minfo.df <- separate(data = minfo.df, col = "info", into = c("mname", "SNP_pos", "align_strand"), sep = ":", remove = T)
head(minfo.df, n = 2)
tail(minfo.df, n = 3)
str(minfo.df)

# Add species name to mname
minfo.df$mname <- paste0(species, "_", minfo.df$mname)


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
seq_and_minfo.df[1:30, c("ref", "alt", "ref_nuc")]
# Note: the 'ref_nuc' is from the genome itself, the ref and alt alleles are defined by the genotyping VCF 
#  and as a result, the 'ref' doesn't always equal 'ref_nuc' 
#  but generally speaking the ref_nuc should be either the ref or alt
#  Custom markers may be different, as they were not genotyped using this ref genome

# Add values into the constant columns
seq_and_minfo.df$strand <- rep("NA", times = nrow(seq_and_minfo.df))
seq_and_minfo.df$mtype <- rep("SNP", times = nrow(seq_and_minfo.df))
seq_and_minfo.df$priority <- rep(2, times = nrow(seq_and_minfo.df))

head(seq_and_minfo.df, n = 1)

# How many didn't match the expectation for either the ref or alt nucl?
table(seq_and_minfo.df$ref_nuc == seq_and_minfo.df$ref | seq_and_minfo.df$ref_nuc == seq_and_minfo.df$alt)
# That should be all TRUE

seq_and_minfo.df[which((seq_and_minfo.df$ref_nuc == seq_and_minfo.df$ref | seq_and_minfo.df$ref_nuc == seq_and_minfo.df$alt)==FALSE), ]

# # Create a vector to indicate the non-matching amplicons
# seq_and_minfo.df$ref_allele_match <- seq_and_minfo.df$ref_nuc == seq_and_minfo.df$ref
# 
# head(x = seq_and_minfo.df, n = 1)


### Confirm no issue with nucleotides and hotspot design
# SPECIAL NOTE: It appears that stacks doesn't always call the genome nucleotide the 'reference allele'
for(i in 1:nrow(seq_and_minfo.df)){
  
  if(seq_and_minfo.df$ref_nuc[i]==seq_and_minfo.df$ref[i]){
    
    seq_and_minfo.df$true_ref[i] <- seq_and_minfo.df$ref[i]
    seq_and_minfo.df$true_alt[i] <- seq_and_minfo.df$alt[i]
    
    # If they aren't correctly oriented to the reference genome, flip them
    
  } else if(seq_and_minfo.df$ref_nuc[i]==seq_and_minfo.df$alt[i]){
    
    seq_and_minfo.df$true_ref[i] <- seq_and_minfo.df$alt[i]
    seq_and_minfo.df$true_alt[i] <- seq_and_minfo.df$ref[i]
    
  } else {
    
    print(paste0("The marker at ", seq_and_minfo.df$mname[i], " is unexpected, adding unknown to this record"))
    
    seq_and_minfo.df$true_ref[i] <- "unkn"
    seq_and_minfo.df$true_alt[i] <- "unkn"
    
  }
}

# How many oddballs failed? 
table(seq_and_minfo.df$true_ref=="unkn")

# # If ok with losing the marker: 
# seq_and_minfo.df <- seq_and_minfo.df[which(seq_and_minfo.df$true_ref!="unkn"), ]

nrow(seq_and_minfo.df)
table(seq_and_minfo.df$true_ref=="unkn")

head(seq_and_minfo.df, n = 1)

# Recombine the seq data back together
seq_and_minfo.df$formatted_seq <- paste0(seq_and_minfo.df$left_seq
                                         , "[", seq_and_minfo.df$true_ref, "/", seq_and_minfo.df$true_alt, "]"
                                         , seq_and_minfo.df$right_seq
                                         )


# Write out results (full) for troubleshooting
write.csv(x = seq_and_minfo.df, file = paste0(output.dir, "/seq_and_minfo_all_data.csv"), quote = F, row.names = F)

seq_and_minfo.df <- seq_and_minfo.df[, c("mname", "chr.x", "pos", "pos"
                               , "true_ref", "true_alt"
                               , "strand", "mtype"
                               , "priority", "formatted_seq"
                               )]

head(seq_and_minfo.df)

write.csv(x = seq_and_minfo.df, file = paste0(output.dir, "/seq_and_minfo_for_submission.csv"), quote = F, row.names = F)



#### Next Steps: #####
# Can do some data checking, and then can submit the markers for designing primers
