# Update run names to remove any instance of double underscores from TVC
# Prepared by Ben J. G. Sutherland (VIU, SBIO)
# Initialized 2026-11-25

# Requires
# - input TVC file

clean_run_ids <- function(tvc_file.FN = NULL){
  
  #### 01. Read in sample data ####
  # Reporting
  print(paste0("Reading in TVC file: ", tvc_file.FN))
  
  # Read in file
  tvc.df <- read.delim2(file = tvc_file.FN, header = T, sep = "\t", check.names = F)
  
  
  #### 02. Clean run names #### 
  ## amplitools requires no double underscores in run names, so convert to single underscores
  
  # Reporting
  print("Updating run names.")
  
  tvc.df$`Run Name` <- gsub(pattern = "__", replacement = "_", x = tvc.df$`Run Name`)
  
  gc()
  
  #### 03. Write out updated TVC files
  ## Set output filenames
  # TVC output
  tvc_updated.FN <- gsub(pattern = "\\.xls", replacement =  "_clean_run_id.xls", x = tvc_file.FN)
  
  # Reporting
  cat("Writing out the updated TVC file as: ", "\n"
      , tvc_updated.FN, "\n"
      )
  
  # Write out updated TVC file
  write.table(x = tvc.df
              , file = tvc_updated.FN
              , row.names = FALSE
              , quote = F
              , sep = "\t"
              )
  
}

