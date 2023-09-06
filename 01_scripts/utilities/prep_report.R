# Prepare report based CKMR output
# Ben J. G. Sutherland
# 2023-08-14

prep_report <- function(relationship = "PO", input.FN = NULL, offspring_ids = NULL){
  
  # If no input.FN is provided, run function on all po*.txt report files in folder
  # note: this approach will likely be removed in the future, only running prep_report from within the parentage wrapper
  if(is.null(input.FN)){
    
    ## Identify files
    files.vec <- NULL
    
    if(relationship=="PO"){
      
      files.vec <- list.files(path = "03_results/", pattern = "^po")
      files.vec <- files.vec[grep(pattern = "\\.txt", x = files.vec)]
      files.vec <- files.vec[grep(pattern = "\\_report.txt$", x = files.vec, invert = T)]
      files.vec <- files.vec[grep(pattern = "\\_no_cutoff.txt$", x = files.vec, invert = T)]
      
    }else if(relationship=="FS"){
      
      print("This output is not yet implemented.")
      
    }
 
  # If an input.FN is provided, and it is a PO result file, use this as the input vector   
  }else if(!is.null(input.FN)){
    
    if(relationship=="PO"){
      
      # Temporary, to match with previous method
      input.FN <- gsub(pattern = "03_results/", replacement = "", x = input.FN)
      
      files.vec <- input.FN
      
    }else{
      
      print("This output is not yet implemented.")
      
    }
  }

  # Reporting
  print("Providing reports for the following result files: ")
  print(paste0("**", files.vec, "**"))
  
  ## Read in files
  # Reporting
  print("Reading in result files")
  
  result.list <- list(); filename <- NULL
  for(i in 1:length(files.vec)){
  
    filename <- files.vec[i]
    
    result.list[[filename]] <- read.delim2(file = paste0("03_results/", files.vec[i]))
    
    # result.list[[filename]][order(result.list[[filename]]$D2_indiv, )]
    
  }
  
  # Debugging
  #print(names(result.list))
  
  ## For each input file, order it, then write out
  for(r in 1:length(result.list)){

    print(paste0("Working on file **", names(result.list[r]), "**"))
    
    # Extract the target data
    data.df <- result.list[[r]]
    
    # Ensure logl ratio is numeric
    data.df$logl_ratio <- as.numeric(data.df$logl_ratio)
    
    # Order by D2, then by decreasing logl ratio
    data.df <- data.df[order(data.df$D2_indiv, data.df$logl_ratio, decreasing = T), ]
    
    # Identify unique offspring
    D2_unique.vec <- unique(x = data.df$D2_indiv)
    
    # If there were offspring IDs provided (i.e., was run as a wrapper)
    if(!is.null(offspring_ids)){ 
      
      # Include any offspring that were completely unassigned
      D2_unique.vec <- unique(c(D2_unique.vec, offspring_ids))
      
      }
    
    # Prepare the output matrix
    output.df <- matrix(data = NA, nrow = length(D2_unique.vec), ncol = 8)
    colnames(output.df) <- c("indiv", "p1", "p1_logl", "p1_nloc", "p2", "p2_logl", "p2_nloc", "other_assigns")
    rownames(output.df) <- D2_unique.vec
    
    output.df <- as.data.frame(output.df)
    
    num_assigns <- NULL; indiv.df <- NULL; 
    for(indiv in 1:nrow(output.df)){
      
      # How many assignments are in the result for this D2 indiv?
      num_assigns <- nrow(data.df[data.df$D2_indiv==D2_unique.vec[indiv],])
      
      # Obtain the names of any extra assignments
      extra_assigns <- NULL; extra_assign_indiv <- NULL; extra_assign_logl <- NULL
      if(num_assigns > 2){
        
        # Obtain elements
        extra_assign_indiv <- data.df[data.df$D2_indiv==D2_unique.vec[indiv], "D1_indiv"][3:num_assigns]
        extra_assign_logl  <- data.df[data.df$D2_indiv==D2_unique.vec[indiv], "logl_ratio"][3:num_assigns]
        extra_assign_logl  <- round(x = as.numeric(extra_assign_logl), digits = 1)
        
        # Build composite
        extra_assigns <- paste0(extra_assign_indiv, "(", extra_assign_logl, ")")
        
      }
      
      # If no assignments at all
      if(num_assigns==0){
        
        # Collect the name of the individual
        output.df[indiv, "indiv"] <- rownames(output.df)[indiv]
        output.df[indiv,"other_assigns"] <- ""
        
      }
      
      # If there are more than zero assignments
      if(num_assigns > 0){
        
        # Obtain the top two assignments
        indiv.df <- head(x = data.df[data.df$D2_indiv==D2_unique.vec[indiv],], n = 2)
        
        output.df[indiv, "indiv"] <- indiv.df[1,"D2_indiv"]
        
        # Obtain info for p1
        output.df[indiv, "p1"]      <- indiv.df[1,"D1_indiv"]
        output.df[indiv, "p1_logl"] <- indiv.df[1,"logl_ratio"]
        output.df[indiv, "p1_nloc"] <- indiv.df[1,"num_loc"]
        
        # Obtain info for p2
        output.df[indiv, "p2"]      <- indiv.df[2,"D1_indiv"]
        output.df[indiv, "p2_logl"] <- indiv.df[2,"logl_ratio"]
        output.df[indiv, "p2_nloc"] <- indiv.df[2,"num_loc"]
        
        # Add any extra assigns
        output.df[indiv, "other_assigns"] <- paste0(extra_assigns, collapse = "; ")

      }
    }
    
    # Convert logls to numeric and round
    output.df$p1_logl <- round(x = as.numeric(output.df$p1_logl), digits = 2)
    output.df$p2_logl <- round(x = as.numeric(output.df$p2_logl), digits = 2)
    
    # Reorder alphanumerically
    output.df <- output.df[order(output.df$indiv), ]
    
    # Write out
    output.FN <- gsub(pattern = "\\.txt$", replacement = "_report.txt", x = names(result.list)[r])
    output.FN <- paste0("03_results/", output.FN)
    print(paste0("Saving output as ", output.FN))
    write.table(x = output.df, file = output.FN, row.names = F, sep = "\t", quote = F)
    
  }
  
}
  