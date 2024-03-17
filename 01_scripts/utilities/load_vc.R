# Load variantcaller data
#  loads all .xls (tab-delim) files in 02_input_data, reduces to only needed cols, 
#   provides number of unique markers, creates "identifier" column, does NOT remove novel variants
#   saves output to input.list
#   "identifier" is comprised of run name, barcode, sample name (all separated by "__")
# Sutherland Bioinformatics, 2023-08-04

load_vc <- function(input_folder = "02_input_data", test_only = FALSE){
  
  inputs <- list.files(path = input_folder, pattern = ".xls")
  
  # Identify sample files for test dataset or for experimental dataset
  if(test_only == FALSE){
    
    # Remove test data from dataset inputs
    inputs <- inputs[inputs!="test_data.xls"]
    
  }else if(test_only == TRUE){
    
    # Reporting
    print(paste0("Test run with ", input_folder, "/test_data.xls"))
    
    # Only use test data
    inputs <- inputs[inputs=="test_data.xls"]
    
  }
  
  # Reporting
  print("Analyzing the following input files: ")
  print(inputs)
  
  # Set nulls
  input.df <- NULL; input.list <- list()
  
  for(file in 1:length(inputs)){
    
    ## Retain filename
    df.name <- gsub(pattern = ".xls", replacement = ".df", x = inputs[file])
    
    ## Read in input file
    # Reporting
    print(paste0("**Reading in file ", inputs[file], "**"))
    
    # Read in
    input.df <- read.delim(file = paste0(input_folder, "/", inputs[file]), header = T, sep = "\t")
    
    # Reporting
    print(paste0("Data loaded as ", nrow(input.df), " rows and ", ncol(input.df), " cols"))
    
    ## Remove extra columns
    print("Reducing columns to retain essentials")
    
    input.df  <- input.df[,c("Chrom", "Position", "Ref"
                             , "Variant", "Allele.Call", "Type", "Allele.Source"
                             , "Allele.Name", "Region.Name", "Coverage", "Strand.Bias"
                             , "Sample.Name", "Barcode", "Run.Name")]
    
    # Reporting
    print(paste0("Data limited to ", nrow(input.df), " rows and ", ncol(input.df), " cols"))
    head(input.df, n = 3)
    print(paste0("Currently, there are ", length(unique(input.df$Allele.Name)), " unique markers"))
    
    
    ## Create identifier
    print("Creating new identifier comprised of Run.Name, Barcode, and Sample.Name, connected by underscores") 
    input.df$identifier <- paste0(input.df$Run.Name, "__", input.df$Barcode, "__", input.df$Sample.Name)
    
    # Reorder cols
    input.df <- input.df %>% select(identifier, everything())
    
    head(input.df, n = 3)
    
    # Reporting
    print("Saving object to input.list")
    input.list[[df.name]] <- input.df
    
    assign(x = "input.list", value = input.list, envir = .GlobalEnv)
  
  }
}
