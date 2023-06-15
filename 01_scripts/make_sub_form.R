# Script to prepare submission form for amplicon panel
# B. Sutherland (Sutherland Bioinformatics)


#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("tidyr")

library("tidyr")

# Set working directory
if(Sys.info()["nodename"] == "Wayne.local"){ 
  print("On Wayne, ready to go")
  setwd("~/Documents/00_sutherland_bioinformatics/GBMF_UBC_Pacific_oyster/ms_oyster_panel/") 
} else if(Sys.info()["nodename"] == "Xavier"){
  print("On Xavier, need to set path")
  #setwd("~/Documents/01_moore_oyster_project/stacks_workflow_all_data/") # Xavier
} else {
  print("You are on an unrecognized system, please set working directory manually")
}


## Info
# sessionInfo()

# Set variables
resources.dir <- "02_input_data//"
input.dir <- "04_extract_loci//"
output.dir <- "05_submission_form//"

species <- "Cgig"

#### 1. Import data ####
### Import sequence data ###
seq.df <- read.delim2(file = paste0(input.dir, "selected_chr_and_seq.txt"), header = F)
seq.df <- as.data.frame(x = seq.df)
colnames(seq.df) <- c("chr_info", "seq")
head(seq.df, n = 1)

# Separate columns into individual details
seq.df <- separate(data = seq.df, col = "chr_info", into = c("mname", "chr_info"), sep = "::", remove = T)
seq.df$mname <- gsub(pattern = "mname_", replacement = "", x = seq.df$mname)
seq.df <- separate(data = seq.df, col = "chr_info", into = c("chr", "pos_range"), sep = ":", remove = T)
seq.df <- separate(data = seq.df, col = "pos_range", into = c("lower_range", "upper_range"), sep = "-", remove = T)

head(x = seq.df, n = 1)

# Update formats
seq.df$seq         <- as.character(seq.df$seq)
seq.df$lower_range <- as.numeric(seq.df$lower_range)
seq.df$upper_range <- as.numeric(seq.df$upper_range)
head(x = seq.df, n = 1)
str(seq.df)


### Import marker data ###
minfo.df <- read.delim2(file = paste0(input.dir, "vcf_selection.csv"), header = F, sep = ",")
minfo.df <- as.data.frame(x = minfo.df, stringsAsFactors = F)
head(minfo.df)
colnames(x = minfo.df) <- c("chr", "pos", "info", "ref", "alt")

minfo.df$ref <- as.character(minfo.df$ref)
minfo.df$alt <- as.character(minfo.df$alt)
minfo.df$chr <- as.character(minfo.df$chr)

head(minfo.df, n = 2)
str(minfo.df)

# Separate columns into individual details
minfo.df <- separate(data = minfo.df, col = "info", into = c("mname", "SNP_pos", "align_strand"), sep = ":", remove = T)
head(minfo.df, n = 2)
str(minfo.df)

#### 02. Combine datasets into one ####
## Merge data
head(x = seq.df, n = 1)
head(minfo.df, n = 1)
nrow(seq.df)
nrow(minfo.df)

seq_and_minfo.df <- merge(x = seq.df, y = minfo.df, by = "mname")
nrow(seq_and_minfo.df)

head(seq_and_minfo.df, n = 2)

# Change nucleotides to all upper case
seq_and_minfo.df$seq <- toupper(x = seq_and_minfo.df$seq)


# Separate the sequence into the first 200 bp, the variant (ref allele), then the second 200 bp 
seq_and_minfo.df$left_seq  <- substr(x = seq_and_minfo.df$seq, start =   1, stop = 200) 
seq_and_minfo.df$ref_nuc   <- substr(x = seq_and_minfo.df$seq, start = 201, stop = 201)
seq_and_minfo.df$right_seq <- substr(x = seq_and_minfo.df$seq, start = 202, stop = 401)

#TODO#: warning, this assumes that all POS are greater than 202 bp from the front, perhaps add a datacheck to verify this #

# Add values into the constant columns
seq_and_minfo.df$strand <- rep("NA", times = nrow(seq_and_minfo.df))
seq_and_minfo.df$mtype <- rep("SNP", times = nrow(seq_and_minfo.df))
seq_and_minfo.df$priority <- rep(2, times = nrow(seq_and_minfo.df))

head(seq_and_minfo.df, n = 1)

# How many didn't match the expectation? 
table(seq_and_minfo.df$ref_nuc == seq_and_minfo.df$ref)

# Create a vector to indicate the non-matching amplicons
seq_and_minfo.df$ref_allele_match <- seq_and_minfo.df$ref_nuc == seq_and_minfo.df$ref

head(x = seq_and_minfo.df, n = 1)


head(seq_and_minfo.df, n = 1)

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

# If ok with losing the marker: 
seq_and_minfo.df <- seq_and_minfo.df[which(seq_and_minfo.df$true_ref!="unkn"), ]

nrow(seq_and_minfo.df)
table(seq_and_minfo.df$true_ref=="unkn")

head(seq_and_minfo.df, n = 1)

# Recombine the seq data back together
seq_and_minfo.df$formatted_seq <- paste0(seq_and_minfo.df$left_seq, "[", seq_and_minfo.df$true_ref, "/", seq_and_minfo.df$true_alt, "]", seq_and_minfo.df$right_seq)


# Write out results (full) for troubleshooting
write.csv(x = seq_and_minfo.df, file = paste0(output.dir, "seq_and_minfo_all_data.csv"), quote = F, row.names = F)

seq_and_minfo.df <- seq_and_minfo.df[, c("mname", "chr.x", "pos", "pos"
                               , "true_ref", "true_alt"
                               , "strand", "mtype"
                               , "priority", "formatted_seq"
                               )]

head(seq_and_minfo.df)

write.csv(x = seq_and_minfo.df, file = paste0(output.dir, "seq_and_minfo_for_submission.csv"), quote = F, row.names = F)



#### Next Steps: #####
# Can do some data checking, and then can submit the markers for designing primers
