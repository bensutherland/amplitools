# Load genepop files
# Sutherland Bioinformatics, 2024-03-17

load_genepops <- function(genepop_folder = "02_input_data/prepped_genepops/", datatype = "SNP"
                          , shorten_name = TRUE){
  
  # Set allele.code
  if(datatype=="SNP"){
    
    allele.code <- 2
    
  }else if(datatype=="microsat"){
    
    allele.code <- 3
    
  }
  
  inputs <- list.files(path = genepop_folder, pattern = ".gen")
  
  
  # Read in all genepop files in input folder, save to my_genepops.list
  my_genepops.list <- list(); shortname <- NULL
  for(i in 1:length(inputs)){
    
    if(shorten_name==TRUE){
      
      # Create shortname for the genepop (first 21 characters)
      shortname <- paste0(substr(inputs[i], start = 1, stop = 21), ".obj")
      
    }else{
      
      shortname <- inputs[i]
      
    }
    
    
    my_genepops.list[[shortname]] <- read.genepop(file = paste0(genepop_folder, inputs[i])
                                                  , ncode = allele.code
                                                  )
    
    print(my_genepops.list[[shortname]])
    
  }
  
  # Reporting
  print("Saved objects to my_genepops.list")
  
  assign(x = "my_genepops.list", value = my_genepops.list, envir = .GlobalEnv)
  
  
}