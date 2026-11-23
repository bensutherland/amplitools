# Update sample names within Torrent VariantCaller file before processing via amplitools
# Prepared by Claire Rycroft (VIU) and Ben J. G. Sutherland (VIU, SBIO)
# Initialized 2026-09-23

# Requires
# - tab-delim sample key with header: original.id \t updated.id

# note: use the following to produce


#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

## Load libraries
library('dplyr')
library('stringr')
library('data.table')


#### FUNCTION FOR BUILDING KEY FILE
df <- fread("data.csv", select = "score")



## Read input data
cat_run1   <- read.delim(file = '02_input_data/G0626-24-VIUN_run1.xls', header = T, sep = "\t")
cat_run2   <- read.delim(file = '02_input_data/G0626-24-VIUN_run2.xls', header = T, sep = "\t")
sample_key <- read.delim(file = '00_archive/pti2026002_sample_id_key.txt')


## Inspect sample names
# Run 1
cat_run1.df <- as.data.frame(cat_run1)
dim(cat_run1.df)
colnames(cat_run1.df)
unique(cat_run1.df$Sample.Name)

# Run 2
cat_run2.df <- as.data.frame(cat_run2)
dim(cat_run2.df)
colnames(cat_run2.df)
unique(cat_run2.df$Sample.Name)


## Use the sample key to replace Sample Names in the TVC file
# sample_id is what is in the original TVC file to be replaced
# true_id is the actual name that the sample should be, corrected by the provider

#table(cat_run1.df$Sample.Name)

cat_run1_orgn.df <- cat_run1.df
cat_run2_orgn.df <- cat_run2.df

# Add a dummy column
# cat_run1.df$Sample.Name.true <- NA
# cat_run2.df$Sample.Name.true <- NA

# # Loop over the rows of the sample key and do a gsub per cat run to replace the incorrect ID with the correct ID
# for(i in 1:nrow(sample_key)){
#   
#   #print(paste0("Converting the string ", sample_key[i,"sample_id"], " to the string ", sample_key[i,"true_id"]))
#   
#   # Fix cat run 1 file
#   cat_run1.df$Sample.Name <- gsub(pattern = sample_key[i,"sample_id"], replacement = sample_key[i,"true_id"], x = cat_run1.df$Sample.Name)
#   
#   # Fix cat run 2 file
#   cat_run2.df$Sample.Name <- gsub(pattern = sample_key[i,"sample_id"], replacement = sample_key[i,"true_id"], x = cat_run2.df$Sample.Name)
#   
# }

# table(cat_run1.df$Sample.Name)

#unique(cat_run1_orgn.df[cat_run1_orgn.df$Sample.Name==cat_run1.df$Sample.Name, "Sample.Name"])

# unique(cat_run1.df$Sample.Name)
# unique(cat_run2.df$Sample.Name)


## Prepare run 2 sample IDs for the rename (requires that there is no preceding text to sample ID)
# Remove string before the target string to use front of string prior to ID as requirement in pattern matching
cat_run2.df$Sample.Name <- gsub(pattern = "G0626-24-VIUN_p1_|G0626-24-VIUN_p2_|G0626-24-VIUN_p3_|G0626-24-VIUN_p4_", replacement = "", x = cat_run2.df$Sample.Name)


## Use the key to replace the old name with the new name
# Loop over the rows of the sample key and do a gsub per cat run to replace the incorrect ID with the correct ID
for(i in 1:nrow(sample_key)){

  # Fix cat run 1 file
  cat_run1.df$Sample.Name <- gsub(pattern = paste0("^", sample_key[i,"sample_id"]), replacement = paste0("new_", sample_key[i,"true_id"]), x = cat_run1.df$Sample.Name
                                  , perl = T)

  # Fix cat run 2 file
  cat_run2.df$Sample.Name <- gsub(pattern = paste0("^", sample_key[i,"sample_id"]), replacement = paste0("new_", sample_key[i,"true_id"]), x = cat_run2.df$Sample.Name
                                  , perl = T)
  

}


## Remove the 'new' prefix (used to identify which samples had been already replaced)
cat_run1.df$Sample.Name <- gsub(pattern = "new_", replacement = "", x = cat_run1.df$Sample.Name)
cat_run2.df$Sample.Name <- gsub(pattern = "new_", replacement = "", x = cat_run2.df$Sample.Name)

# Check to see the original and new name for those that changed during the rename
unique(
  paste0(cat_run1.df$Sample.Name[cat_run1.df$Sample.Name!=cat_run1_orgn.df$Sample.Name]
         , "__"
         , cat_run1_orgn.df$Sample.Name[cat_run1.df$Sample.Name!=cat_run1_orgn.df$Sample.Name])
  )

### DID NOT WORK ###
# cat_run1_orgn.df$Sample.Name[grep(pattern = sample_key[i,"sample_id"], x = cat_run1_orgn.df$Sample.Name) & !grep(pattern = "new_", x = cat_run1_orgn.df$Sample.Name)]
# 
# result <- text_vector[grepl("a", text_vector) & !grepl("r", text_vector)]
### /END/ DID NOT WORK ###


### ORIGINAL ###

# # Update names from key
# ## CAT run 1
# cat_run1_corrected <- cat_run1 %>%
#   left_join(.
#             , sample_key
#             , by = c('Sample.Name' = 'sample_id')
#             ) %>%
#   # confirm changes
#   # mutate(Sample.Name_corrected = case_when(Sample.Name != true_id ~ true_id
#   #                                          , TRUE ~ Sample.Name)) %>%
#   # select(Sample.Name
#   #        , Sample.Name_corrected
#   #        , true_id
#   #        , everything())
#   mutate(Sample.Name = case_when(Sample.Name != true_id ~ true_id
#                                  , TRUE ~ Sample.Name)) %>%
#   select(!true_id)
# 

# ## CAT run 2
# cat_run2_corrected <- cat_run2
# 
# for (i in seq_len(nrow(sample_key))) {
#   
#   edit <- str_detect(cat_run2_corrected$Sample.Name
#                      , fixed(sample_key$sample_id[i]))
#   
#   cat_run2_corrected$Sample.Name_corrected[edit] <- str_replace(cat_run2_corrected$Sample.Name[edit]
#                                                           , fixed(sample_key$sample_id[i])
#                                                           , sample_key$true_id[i]
#                                                           )
# }
# 
# ## confirm dataframe changes
# # cat_run2_corrected <- cat_run2_corrected %>% 
# #   select(Sample.Name
# #          , Sample.Name_corrected
# #          , everything())
# 
# ## Restore column names
# cat_run2_corrected <- cat_run2_corrected %>% 
#   select(!Sample.Name) %>% 
#   rename(Sample.Name = Sample.Name_corrected)

#### Write output ####
# Write corrected dataframes to file
write.table(cat_run1.df
            , '03_results/G0626-24-VIUN_run1_corrected.xls'
            , row.names = FALSE
            , quote = F
            , sep = "\t"
            )

write.table(cat_run2.df
            , '03_results/G0626-24-VIUN_run2_corrected.xls'
            , row.names = FALSE
            , quote = F
            , sep = "\t"
            )

# end #
