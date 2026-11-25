# Update sample names within Torrent VariantCaller (TVC) file before processing via amplitools
# Prepared by Claire Rycroft (VIU) and Ben J. G. Sutherland (VIU, SBIO)
# Initialized 2026-09-23

# Requires
# - tab-delim sample key with header: original.id \t updated.id

correct_tvc_ids <- function(sample_key.FN = NULL){
  
  #### 01. Read in sample data ####
  # Reporting
  print(paste0("Reading in: ", sample_key.FN))
  
  ## Read in sample key
  sample_key.df <- read.delim2(file = sample_key.FN, header = T, sep = "\t")
  head(sample_key.df)
  
  ## Add replicate ID for duplicates
  # Check to see if there are any instances of more than one duplicate ID
  if(any(table(sample_key.df$updated.id) > 2)){
    
    stop("This program has not implemented handling more than two duplicate values in updating names.")
    
  }else {
    
    print("Duplicate target IDs detected, adding '_rep' to second instance.")
    
  }
  
  # Add the replicate indicator for second instance of identifier
  sample_key.df[duplicated(sample_key.df$updated.id), "updated.id"] <- paste0(sample_key.df[duplicated(sample_key.df$updated.id), "updated.id"], "_rep")
  
  
  ## Read in sample file
  # Automatically set filename
  tvc_raw.FN <- gsub(pattern = "_sample_key_annot.txt", replacement = ".xls", x = sample_key.FN)
  
  # Reporting
  print(paste0("Reading in TVC file: ", tvc_raw.FN))
  
  # Read in file
  tvc_raw.df <- read.delim2(file = tvc_raw.FN, header = T, sep = "\t")
  
  
  #### 02. Replace sample names #### 
  ## Use the sample key to replace the old name with the new name
  
  # Reporting
  print("Updating existing sample names to target names.")
  
  # Loop over the rows of the sample key and do a gsub per cat run to replace the incorrect ID with the correct ID
  for(i in 1:nrow(sample_key.df)){
  
    # Fix cat run 1 file
    tvc_raw.df$Sample.Name <- gsub(pattern = paste0("^", sample_key.df[i, "original.id"])
                                   , replacement = paste0("new_", sample_key.df[i,"updated.id"])
                                   , x = tvc_raw.df$Sample.Name
                                   , perl = T
                                   )
  
  
  
  }
  
  ## Remove the 'new' prefix (used to identify which samples had been already replaced)
  tvc_raw.df$Sample.Name <- gsub(pattern = "new_", replacement = "", x = tvc_raw.df$Sample.Name)
  tvc_updated.df <- tvc_raw.df
  
  # Clean up
  rm(tvc_raw.df)
  gc()
  
  #### 03. Write out updated TVC files
  ## Set output filenames
  # TVC output
  tvc_updated.FN <- gsub(pattern = "tvc_to_update\\/", replacement =  "", x = tvc_raw.FN)
  tvc_updated.FN <- gsub(pattern = "\\.xls", replacement =  "_updated_ids.xls", x = tvc_updated.FN)
  
  # Applied sample key
  sample_key_updated.FN <- gsub(pattern = "tvc_to_update\\/", replacement = "", x = sample_key.FN)
  sample_key_updated.FN <- gsub(pattern = "\\.txt", replacement = "_applied.txt", x = sample_key_updated.FN)
  
  # Reporting
  cat("Writing out the updated TVC file as: ", "\n"
      , tvc_updated.FN, "\n"
      , "...and updated sample key as: ", "\n"
      , sample_key_updated.FN
      )
  
  # Write out updated TVC file
  write.table(x = tvc_updated.df
              , file = tvc_updated.FN
              , row.names = FALSE
              , quote = F
              , sep = "\t"
              )
  
  # Write out applied sample key
  write.table(x = sample_key.df
              , file = sample_key_updated.FN
              , row.names = F
              , col.names = T
              , quote = F
              , sep = "\t"
  )

}

