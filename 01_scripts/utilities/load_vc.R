# Load variantcaller data
# Sutherland Bioinformatics, 2023-08-04

load_vc <- function(input_folder = "02_input_data"){
  
  inputs <- list.files(path = input_folder, pattern = ".xls")
  
  # Set nulls
  input.df <- NULL
  
  for(file in 1:length(inputs)){
    
    ## Retain filename
    df.name <- gsub(pattern = ".xls", replacement = ".df", x = inputs[file])
    
    ## Read in input file
    # Reporting
    print(paste0("Reading in file ", inputs[file]))
    
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
    print(paste0("Saving object as ", df.name))
    
    assign(x = df.name, value = input.df, envir = .GlobalEnv)
  
  }
}

  