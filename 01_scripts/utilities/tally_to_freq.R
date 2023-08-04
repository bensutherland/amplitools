# General inspection of variant frequency
# Sutherland Bioinformatics, 2023-08-04
# caution: these are all approximations of MAF, since non-target reference alleles are not genotyped by the software

tally_to_freq <- function(df = "df", allele_source = "novel"){
  
  # Reporting
  print(paste0("Assessing ", allele_source, " variants"))
  
  # Set target of function
  if(allele_source=="novel"){
    
    allele_source <- "Novel"
    type          <- "Frequency of variant"
    
  }else if(allele_source=="hotspot"){
    
    allele_source <- "Hotspot"
    type          <- "Minor allele frequency"
    
  }else(
    
    stop("The allele_source must be either novel or hotspot")
    
  )
  
  # Provide warning
  cat("Caution: \nfrequencies only approximated, \nassumes no missing data and \nabsence of calls as homo ref for novel variants")
  
  # Reporting
  # print("Working with the following file")
  # print(head(df, n = 3))
   
  # Assess how many unique individuals are present in the df
  num_indiv <- length(unique(x = df$identifier))
  print(paste0("There are a total of ", num_indiv, " individuals present in the data"))
  
  # Reporting
  print("Assigning universal locus names (i.e., locus.id) to novel variants for broader comparison")
  
  # Assign locus names to all variants
  df$locus.id <- paste0(df$Chrom, "__", df$Position)
  #print(head(df, n = 3))
  
  # Reporting
  print(paste0("There are a total of ", length(unique(df$locus.id)), " unique loci"))
  
  # Only keep specified variants
  print(tolower(paste0("Retaining only specified variants from the ", allele_source, " allele source")))
  df <- df[df$Allele.Source==allele_source,]
  print(paste0("Retained  ", length(unique(df$locus.id)), " unique loci"))
  
  #print(head(df, n = 3))
  
  # Per row, tally minor allele variants
  print("Tallying variants")
  df$allele.count <- NA
  df[df$Allele.Call=="Homozygous", "allele.count"] <- 2
  df[df$Allele.Call=="Heterozygous", "allele.count"] <- 1
  df[df$Allele.Call=="Absent"|df$Allele.Call=="No Call", "allele.count"] <- 0
  
  #print(head(df, n = 3))
  
  # Tally minor allele (hotspot) or non-ref (novel) variants
  tallies <- df %>% group_by(locus.id) %>%
                  summarise(total_views = sum(allele.count))
  
  tallies <- as.data.frame(tallies)
  
  tallies <- tallies[order(tallies$total_views, decreasing = T),]
  
  #head(tallies)
  
  # Calculate frequency of , assuming 100% genotyping (no missing data)
  tallies$freq <- tallies$total_views / (num_indiv*2)
  
  #head(tallies)
  
  # Convert to minor allele freq
  for(i in 1:nrow(tallies)){
    
    if(tallies[i,"freq"] > 0.5){
      
      tallies[i, "freq"] <- 1 - tallies[i, "freq"]
      
    }
  }
  
  # Reorder again
  tallies <- tallies[order(tallies$total_views, decreasing = T),]
  
  # Determine how many loci are greater than 0.01
  nloc_greater_than_0.01 <- nrow(tallies[tallies$freq > 0.01, ])
  nloc_greater_than_0.1 <- nrow(tallies[tallies$freq > 0.1, ])
  print(paste0("Number loci with MAF > 0.1 = ", nloc_greater_than_0.1))
  
  output_figure.FN <- tolower(paste0("03_results/freq_of_MAF_for_obs_minor_allele_", allele_source, ".pdf"))
  output_df.FN     <- tolower(paste0("03_results/tallies_and_MAF_", allele_source, ".txt"))
  
  
  print(paste0("Plotting output figure as ", output_figure.FN))
  pdf(file = output_figure.FN, width = 6, height = 4)
  hist(tallies$freq, las = 1, breaks = 20, main = ""
       , xlab = "MAF (observed minor allele)"
       )
  
  text(x = 0.3, y = (nrow(tallies) / 3.2), labels = paste0("n = ", nrow(tallies), " loci"))
  text(x = 0.3, y = (nrow(tallies) / 3.5),   labels = paste0("MAF > 0.01 = " , nloc_greater_than_0.01, " loci"))
  dev.off()
  
  print(paste0("Writing out data as ", output_df.FN))
  write_delim(x = tallies, file = output_df.FN, delim = "\t")
  
  # TODO: compare directly MAF for alternates vs. target using the region name (SP)
  
}
