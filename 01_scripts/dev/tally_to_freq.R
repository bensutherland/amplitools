# General inspection of variant frequency
# Sutherland Bioinformatics, 2023-08-04

# Note: the development of this script was halted because one cannot assume that the novel and the hotspot variants within the same Region.Name are genotyped in the same number of individuals, as the genotyping varies across the region. Therefore, one cannot compare frequencies between these types of variants, even within a region. 

# There are currently no plans to continue development of this script. 


tally_to_freq <- function(df = "df", allele_source = "novel"){
  
  # # Match software default capitalization
  # allele_source <- gsub(pattern = "novel", replacement = "Novel", x = allele_source)
  # allele_source <- gsub(pattern = "hotspot", replacement = "Hotspot", x = allele_source)
  # 
  # # Reporting
  # print(paste0("Assessing ", allele_source, " variants"))
  
  # Provide warning
  cat("Caution: \nnovel variant frequencies only approx. (absence of call is assumed reference) \nand AF is calculated based on amplicon's hotspot presence")
  head(df, n = 3)
  
  ## Assign locus names to all loci
  print("Assigning universal locus names (i.e., locus.id) to novel variants for broader comparison")
  df$locus.id <- paste0(df$Chrom, "__", df$Position)
  head(df, n = 3)
  
  ## How many unique individuals are present in the df?
  print(paste0("There are a total of ", length(unique(x = df$identifier)), " individuals present in the data"))
  
  ## How many unique loci are present in the df?
  print(paste0("There are a total of ", length(unique(df$locus.id)), " unique loci"))
  
  # # Retain backup before filtering
  # df.bck <- df

  
  #### 01. Count number of indiv genotyped per amplicon (use hotspot) ####
  # Retain hotspot rows only
  df_hotspot <- df[df$Allele.Source=="Hotspot", ]
  
  # Reporting
  print(paste0("Number of unique hotspot loci: ", length(unique(df_hotspot$locus.id))))
  print(paste0("Number of unique amplicon regions: ", length(unique(df_hotspot$Region.Name))))
  print("Note: if the number of hotspot loci > number of regions, there may be multiple hotspots per region")
  
  # Remove any hotspots that are fully missing data
  print("Removing fully untyped hotspots (i.e., all 'No Call')")
  df_hotspot <- df_hotspot[df_hotspot$Allele.Call!="No Call", ] 
  print(paste0("Number of unique hotspot loci: ", length(unique(df_hotspot$locus.id))))
  
  # Create a vector of locus.ids for genotyped amplicons (must be genotyped in at least one individual)
  total_locus.id <- unique(df_hotspot$locus.id)
  print(paste0("Inspecting ", length(total_locus.id), " amplicon hotspots"))
  
  # Per locus.id, determine how many samples were genotyped
  locus_of_interest <- NULL; region_name <- NULL; num_indiv_genod_at_locus <- NULL; num_indiv_genod_at_locus.vec <- NULL
  for(i in 1:length(total_locus.id)){
    
    locus_of_interest <- total_locus.id[i]
    region_name <- unique(df_hotspot[df_hotspot$locus.id==locus_of_interest, "Region.Name"])
    
    num_indiv_genod_at_locus     <- length(df_hotspot[df_hotspot$locus.id==locus_of_interest, "Sample.Name"])
    num_indiv_genod_at_locus     <- paste0(locus_of_interest, "___", region_name, "___", num_indiv_genod_at_locus)
    num_indiv_genod_at_locus.vec <- c(num_indiv_genod_at_locus, num_indiv_genod_at_locus.vec)
    
  }
  
  # Make into df
  num_indiv_genod_at_locus.df  <- as.data.frame(num_indiv_genod_at_locus.vec)
  head(num_indiv_genod_at_locus.df)
  
  
  # Separate back into locus.id, region.name and num indiv geno'd at locus
  num_indiv_genod_at_locus.df  <- separate(data = num_indiv_genod_at_locus.df
                                           , col = "num_indiv_genod_at_locus.vec"
                                           , into = c("locus.id", "Region.Name", "indiv.genod")
                                           , sep = "___"
                                           , remove = T
                                          )
  
  head(num_indiv_genod_at_locus.df)
  tail(num_indiv_genod_at_locus.df)
  nrow(num_indiv_genod_at_locus.df) # 587 rows
  
  # Ensure numeric
  num_indiv_genod_at_locus.df$indiv.genod <- as.numeric(num_indiv_genod_at_locus.df$indiv.genod)
  
  # Calculate the total genotyped alleles per locus from the total number of indiv x 2 (assumes diploid)
  num_indiv_genod_at_locus.df$total.poss.alleles <- num_indiv_genod_at_locus.df$indiv.genod * 2
  head(num_indiv_genod_at_locus.df)
  
  # In case of multiple hotspots per region, keep the record with the highest geno rate
  num_indiv_genod_at_locus.df     <-  num_indiv_genod_at_locus.df[order(num_indiv_genod_at_locus.df$Region.Name
                                                                        , num_indiv_genod_at_locus.df$indiv.genod
                                                                        , decreasing = T), ]
  
  num_indiv_genod_at_locus.df <- num_indiv_genod_at_locus.df[!duplicated(x = num_indiv_genod_at_locus.df$Region.Name), ]
  
  head(num_indiv_genod_at_locus.df)
  nrow(num_indiv_genod_at_locus.df) # 586 rows
  min(num_indiv_genod_at_locus.df$indiv.genod) # minimum = 1
  max(num_indiv_genod_at_locus.df$indiv.genod) # maximum = 191
  
  ## This df contains number indiv genod per Region.Name
  
  
  #### ISSUE #####
  # 
  # # Bring this information back into the df object, which contains novel and hotspot
  # df_w_num_genod.df  <- merge(x = df, y = num_indiv_genod_at_locus.df, by = "Region.Name", all.x = T)
  # # Rename the df locus.id column (locus.id.x) as locus.id, as the locus.id.y is the hotspot
  # colnames(df_w_num_genod.df)[colnames(df_w_num_genod.df)=="locus.id.x"] <- "locus.id"
  # 
  # print(paste0("Number of rows in input df = ", nrow(df)))
  # print(paste0("Number of rows in df with new info added = ", nrow(df_w_num_genod.df)))
  # 
  # # Rename as df
  # if(nrow(df)==nrow(df_w_num_genod.df)){
  #   
  #   df <- df_w_num_genod.df
  #   
  # }else{
  #   
  #   stop("There are unequal records before and after bringing in the number indiv per locus info")
  #   
  # }
  # 
  # head(df)
  
  
  #### 02. Tally variants ####
  # # Only keep specified variants
  # print(paste0("Retaining only specified variants from the ", tolower(allele_source), " allele source"))
  # df <- df[df$Allele.Source==allele_source,]
  # print(paste0("Retained ", length(unique(df$locus.id)), " unique loci"))
  # 
  # head(df, n = 3)
  
  # Per row, tally minor allele variants
  head(df, n = 3)
  print("Tallying variants")
  df$allele.count <- NA
  df[df$Allele.Call=="Homozygous", "allele.count"] <- 2
  df[df$Allele.Call=="Heterozygous", "allele.count"] <- 1
  df[df$Allele.Call=="Absent"|df$Allele.Call=="No Call", "allele.count"] <- 0
  
  head(df, n = 5) 
  
  #table(df$Allele.Call)
  
  # Retain the Region.Name with the locus ID when tallying
  df$locus.id___Region.Name <- paste0(df$locus.id, "___", df$Region.Name)
  head(df)
  
  # Tally minor allele occurrences by locus ID
  # tallies <- df %>% group_by(locus.id) %>%
  #                 summarise(total_views = sum(allele.count))
  
  # Tally minor allele occurrences by locus ID
  tallies <- df %>% group_by(locus.id___Region.Name) %>%
                summarise(total_views = sum(allele.count))
   
  tallies <- as.data.frame(tallies)
  
  tallies <- separate(data = tallies, col = "locus.id___Region.Name", into = c("locus.id", "Region.Name"), sep = "___", remove = T)
  
  # tallies <- tallies[order(tallies$total_views, decreasing = T),]
  head(tallies)
  str(tallies)
  nrow(tallies)
  head(num_indiv_genod_at_locus.df)
  str(num_indiv_genod_at_locus.df)
  nrow(num_indiv_genod_at_locus.df)

  
  #### MERGE TALLIES WITH NUM INDIV GENOD AT LOCUS
  # nrow(tallies)
  # nrow(num_indiv_genod_at_locus.df)
  all_data.df <- merge(x = tallies, y = num_indiv_genod_at_locus.df, by = "Region.Name", all.x = T)
  nrow(all_data.df) # 11,046
  
  head(all_data.df)
  
  all_data.df[all_data.df$total_views > all_data.df$total.poss.alleles, ]
  
  # Calculate MAF
  all_data.df$MAF <- all_data.df$total_views / all_data.df$total.poss.alleles
  hist(all_data.df$MAF)
  write_delim(x = df, file = "03_results/output.txt", delim = "\t")
  
  #### ISSUE HERE ####
  
  # # Add tallies to df
  # head(df)
  # head(tallies)
  # 
  # nrow(df)
  # all_data.df <- merge(x = df, y = tallies, by = "locus.id")
  # nrow(all_data.df)
  # 
  # if(nrow(df)==nrow(all_data.df)){
  #   
  #   df <- all_data.df
  #   
  # }else{
  #   
  #   stop("The input df and the output df are not of equal size, quitting")
  #   
  # }
  # 
  
  # Combine num_indiv_genod_at_locus.df and tallies based on locus.id
  head(num_indiv_genod_at_locus.df)
  head(tallies)
  length(intersect(x = num_indiv_genod_at_locus.df$locus.id, y = tallies$locus.id)) # 586
  nrow(num_indiv_genod_at_locus.df)
  nrow(tallies)
  
  
  
  # Calculate MAF
  df$MAF <- df$total_views / df$total.poss.alleles
  
  hist(df$MAF)
  
  write_delim(x = df, file = "03_results/output.txt", delim = "\t")
  
  
  # # Convert to minor allele freq
  # for(i in 1:nrow(tallies)){
  #   
  #   if(tallies[i,"freq"] > 0.5){
  #     
  #     tallies[i, "freq"] <- 1 - tallies[i, "freq"]
  #     
  #   }
  # }
  # 
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
