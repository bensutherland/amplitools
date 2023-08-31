# Visualize network of relatedness 
# Sutherland Bioinformatics
# 2023-08-02

graph_relatives <- function(input.FN = "03_results/po_broodstock_vs_spat_pw_logl_5.txt"
                            , logl_cutoff = 5
                            , drop_string = "G00" 
                            , directed = FALSE
                            , plot_width = 5
                            , plot_height = 5
                            ){
  
  # Read in input
  interacts.df <- read.delim2(file = input.FN
                        , header = T, sep = "\t"
  )
  
  # Format the logl ratio as numeric
  interacts.df$logl_ratio <- as.numeric(interacts.df$logl_ratio)
  
  # Reduce length of names if drop_string is set
  interacts.df$D2_indiv <- gsub(pattern = drop_string, replacement = "", x = interacts.df$D2_indiv)
  interacts.df$D1_indiv <- gsub(pattern = drop_string, replacement = "", x = interacts.df$D1_indiv)
  
  # Reporting
  print(paste0("There are a total of ", nrow(interacts.df), " pairs in the input")) # how many interactions? 
  
  # Reporting
  print(paste0("Reducing pairs to only those with log likelihood value > ", logl_cutoff))
  
  # Reduce to only those greater than the logl
  interacts.df <- interacts.df[interacts.df$logl_ratio > logl_cutoff, ]
  
  # Reporting
  print(paste0("After filtering pairs, there are a total of ", nrow(interacts.df), " pairs")) # how many interactions? 
  
  # Create plot
  interacts.plot <- graph_from_data_frame(interacts.df, directed = directed)
  #print(interacts.plot, e=TRUE, v=TRUE)
  
  # plot.igraph(interacts.plot
  #             , label.cex = 0.2
  #             , vertex.color = "grey"
  #             , vertex.size = 10
  #             , vertex.label.color = "black"
  #             , vertex.frame.color="gray"
  #             , vertex.label.cex=0.8
  #             , vertex.label.dist=2
  #             , edge.curved=0.2
  # )
  
  # Set output filename
  out_plot.FN <- gsub(pattern = ".txt", replacement = paste0("_filtered_to_logl_", logl_cutoff, "_plot.pdf"), x = input.FN)
  
  # Plot and output
  pdf(file = out_plot.FN, width = plot_width, height = plot_height)
  plot.igraph(interacts.plot
              , vertex.color = "grey"
              , vertex.label.color = "black"
              , vertex.frame.color="gray"
              , vertex.size = 5
              , label.cex = 1
              , vertex.label.cex=0.7
              , vertex.label.dist=0.5
              , edge.curved=0.2
              , edge.color = "black"
              , edge.width = 1
  )
  dev.off()

}






