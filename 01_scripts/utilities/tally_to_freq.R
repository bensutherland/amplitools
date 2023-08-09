# General inspection of variant frequency
# Sutherland Bioinformatics, 2023-08-04
# caution: these are all approximations of MAF, since non-target reference alleles are not genotyped by the software

tally_to_freq <- function(df = "df", allele_source = "novel"){
  
  # Replace lower case with the  capitalized term for later matching
  allele_source <- gsub(pattern = "novel", replacement = "Novel", x = allele_source)
  allele_source <- gsub(pattern = "hotspot", replacement = "Hotspot", x = allele_source)
  
  # Reporting
  print(paste0("Assessing ", allele_source, " variants"))
  
  # Provide warning
  cat("Caution: \nnovel variant frequencies only approx. (absence of call is assumed reference) \nand AF is calculated based on amplicon's hotspot presence")
  head(df, n = 3)
  
  ## Assign locus names to all variants
  print("Assigning universal locus names (i.e., locus.id) to novel variants for broader comparison")
  df$locus.id <- paste0(df$Chrom, "__", df$Position)
  head(df, n = 3)
  
  ## How many unique individuals are present in the df?
  print(paste0("There are a total of ", length(unique(x = df$identifier)), " individuals present in the data"))
  
  ## How many unique loci are present in the df?
  print(paste0("There are a total of ", length(unique(df$locus.id)), " unique loci"))
  
  # Retain backup before filtering
  df.bck <- df

  #### 01. Num indiv genotyped at each amplicon ####
  ## Determine how many individuals are genotyped at each region (i.e., amplicon)
  # Retain only hotspot info
  df_hotspot <- df[df$Allele.Source=="Hotspot", ]
  head(df_hotspot)
  nrow(df_hotspot)
  # Remove any sample-locus that is a no call (i.e., missing)
  df_hotspot <- df_hotspot[df_hotspot$Allele.Call!="No Call", ] # drop no calls
  nrow(df_hotspot)
  
  # What are the unique locus IDs? Create a vector
  total_locus.id <- unique(df_hotspot$locus.id)
  length(total_locus.id)
  
  # Loop over the vector, and for each unique locus ID, determine how many samples were typed at that locus
  # note: this may be replaced by summarise (later)
  locus_of_interest <- NULL; num_indiv_genod_at_locus <- NULL; num_indiv_genod_at_locus.vec <- NULL
  for(i in 1:length(total_locus.id)){
    
    locus_of_interest <- total_locus.id[i]
    
    num_indiv_genod_at_locus     <- length(df_hotspot[df_hotspot$locus.id==locus_of_interest, "Sample.Name"])
    num_indiv_genod_at_locus     <- paste0(locus_of_interest, "___", num_indiv_genod_at_locus)
    num_indiv_genod_at_locus.vec <- c(num_indiv_genod_at_locus.vec, num_indiv_genod_at_locus)
    
  }
  
  head(num_indiv_genod_at_locus.vec)
  
  # Make into df
  num_indiv_genod_at_locus.df  <- as.data.frame(num_indiv_genod_at_locus.vec)
  num_indiv_genod_at_locus.df  <- separate(data = num_indiv_genod_at_locus.df, col = "num_indiv_genod_at_locus.vec", into = c("locus.id", "indiv.genod")
                                                    , sep = "___", remove = T
                                          )
  
  num_indiv_genod_at_locus.df$indiv.genod <- as.numeric(num_indiv_genod_at_locus.df$indiv.genod)
  num_indiv_genod_at_locus.df$total.poss.alleles <- num_indiv_genod_at_locus.df$indiv.genod * 2
  head(num_indiv_genod_at_locus.df)
  
  merge(x = num_indiv_genod_at_locus.df, y = df, by = "locus.id")
  
  
  #### 02. Tally variants ####
  # Only keep specified variants
  print(paste0("Retaining only specified variants from the ", tolower(allele_source), " allele source"))
  df <- df[df$Allele.Source==allele_source,]
  print(paste0("Retained ", length(unique(df$locus.id)), " unique loci"))
  
  head(df, n = 3)
  
  # Per row, tally minor allele variants
  print("Tallying variants")
  df$allele.count <- NA
  df[df$Allele.Call=="Homozygous", "allele.count"] <- 2
  df[df$Allele.Call=="Heterozygous", "allele.count"] <- 1
  df[df$Allele.Call=="Absent"|df$Allele.Call=="No Call", "allele.count"] <- 0
  
  head(df, n = 3) 
  
  table(df$Allele.Call)
  
  # Tally minor allele (hotspot) or non-ref (novel) variants
  tallies <- df %>% group_by(locus.id) %>%
                  summarise(total_views = sum(allele.count))
  
  tallies <- as.data.frame(tallies)
  
  tallies <- tallies[order(tallies$total_views, decreasing = T),]
  
  head(tallies)
  
  
  #### Bring it all together ####
  head(num_indiv_genod_at_locus.df) # this is actually the hotspot itself
  head(tallies)
  
  
  
  # Obtain the Region.Name for each locus.id to attach to tallies
  locus_and_regions <- df.bck[,c("locus.id", "Region.Name")]
  head(locus_and_regions)
  locus_and_regions <- paste0(locus_and_regions$locus.id, "___", locus_and_regions$Region.Name)
  locus_and_regions <- unique(locus_and_regions)
  locus_and_regions <- as.data.frame(locus_and_regions)
  head(locus_and_regions)
  locus_and_regions <- separate(data = locus_and_regions, col = "locus_and_regions", into = c("locus", "region"), sep = "___", remove = T)
  head(locus_and_regions)
  
  test <- merge(x = tallies, y = locus_and_regions, by.x = "locus.id", by.y = "locus", all.x = T)
  head(test)
  nrow(test)
  
  # Now have the region, the locus id, and the total number of views
  
  
  
  # Combine with tallies
  head(tallies)
  nrow(tallies)
  
  all_data.df <- merge(x = tallies, y = num_indiv_genod_at_locus.df, by.x = "locus.id", by.y = "locus.id")
  head(all_data.df)
  nrow(all_data.df)
  
  #plot(all_data.df$total_views, all_data.df$total.poss.alleles)
  
  all_data.df$total_views / all_data.df$total.poss.alleles
  
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
