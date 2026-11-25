# Prepare a sample key to use with the TVC file samplename updater script
# Prepared by Ben J. G. Sutherland (VIU, SBIO)
# Initialized 2026-11-23

# Requires
# - TVC files to be updated in 02_input_data/tvc_to_update/

prep_sample_key <- function(){

  # Set path to raw TVC file(s) to be updated (do not change)
  input_raw_TVC.path <- "02_input_data/tvc_to_update"
  
  ## Prepare sample key file(s), one for each TVC file in source folder
  # Identify TVC files in the path
  raw_tvc_files.vec <- list.files(path = input_raw_TVC.path, pattern = ".xls", full.names = T)
  
  # For each raw TVC file, prepare a blank sample key file with all sample names shown to be filled
  sample_ids.vec <- NULL; tvc_blank_key.df <- NULL 
  for(i in 1:length(raw_tvc_files.vec)){
    
    # Read in file 
    sample_ids.vec <- fread(file = raw_tvc_files.vec[i], sep = "\t", header = T, select = "Sample Name")
    
    # Pull sample IDs into vector
    sample_ids.vec <- unique(sample_ids.vec$`Sample Name`)
    
    # Convert to df
    tvc_blank_key.df$original.id <- sample_ids.vec
    
    # Add empty column for filling
    tvc_blank_key.df$updated.id <- NA
    
    # Prep output filename
    output.FN <- gsub(pattern = "\\.xls$", replacement = "_sample_key.txt", x = raw_tvc_files.vec[i])
    
    # Write out
    write.table(x = tvc_blank_key.df, file = output.FN, sep = "\t"
                , row.names = F, col.names = T, quote = F
    )
    
    }

  # Reporting 
  print(paste0("The following sample keys have been written to 02_input_data/00_tvc_raw: " ))
  cat(raw_tvc_files.vec, "\n")
  
  # Reporting
  print("Next: fill in the updated.id column in each *_sample_key.txt and save as *_sample_key_annot.txt in the same folder")
  
}
