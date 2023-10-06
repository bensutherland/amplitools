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
  df_hotspot <- df_hotspot[df_hotspot$Allele.Call!="No Call", ] # Remove rows where the allele call is no call
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
  nrow(num_indiv_genod_at_locus.df)
  
  # Ensure numeric
  num_indiv_genod_at_locus.df$indiv.genod <- as.numeric(num_indiv_genod_at_locus.df$indiv.genod)
  
  # Derive number alleles (total number of indiv x 2; assumes diploid)
  num_indiv_genod_at_locus.df$total.poss.alleles <- num_indiv_genod_at_locus.df$indiv.genod * 2
  head(num_indiv_genod_at_locus.df)
  
  # In case of multiple hotspots per region, keep the record with the highest geno rate
  num_indiv_genod_at_locus.df     <-  num_indiv_genod_at_locus.df[order(num_indiv_genod_at_locus.df$Region.Name
                                                                        , num_indiv_genod_at_locus.df$indiv.genod
                                                                        , decreasing = T), ]
  
  num_indiv_genod_at_locus.df <- num_indiv_genod_at_locus.df[!duplicated(x = num_indiv_genod_at_locus.df$Region.Name), ]
  
  # Inspect
  head(num_indiv_genod_at_locus.df)
  nrow(num_indiv_genod_at_locus.df)
  min(num_indiv_genod_at_locus.df$indiv.genod) # What is the minimum number of genotyped samples in the data? 
  max(num_indiv_genod_at_locus.df$indiv.genod) # What is the maximum number of genotyped samples in the data? 
  
  
  #### 02. Connect approx. number geno'd samples per region back to the original df w/ non-hotspot ####
  # Bring this information back into the df object, which contains novel and hotspot
  dim(df) # 369,575 rows, 16 cols
  df_w_num_genod.df  <- merge(x = df, y = num_indiv_genod_at_locus.df, by = "Region.Name", all.x = T)
  dim(df_w_num_genod.df)
  # note: locus.id.y is the type locus, or the hotspot locus for the region
  # note: locus.id.x is the exact locus ID, where it is chr and pos
  
  # Rename the df locus.id column (locus.id.x) as locus.id, as the locus.id.y is the hotspot
  colnames(df_w_num_genod.df)[colnames(df_w_num_genod.df)=="locus.id.x"] <- "locus.id"
  colnames(df_w_num_genod.df)[colnames(df_w_num_genod.df)=="locus.id.y"] <- "locus.id.rel.hotspot"

  print(paste0("Number of rows in input df = ", nrow(df)))
  print(paste0("Number of rows in df with new info added = ", nrow(df_w_num_genod.df)))

  # Rename as df
  if(nrow(df)==nrow(df_w_num_genod.df)){

    df <- df_w_num_genod.df

  }else{

    stop("There are unequal records before and after bringing in the number indiv per locus info")

  }

  head(df)

  
  #### 02. Tally variants ####
  
  ## NOT NEEDED? ##
  # # Only keep specified variants
  # print(paste0("Retaining only specified variants from the ", tolower(allele_source), " allele source"))
  # df <- df[df$Allele.Source==allele_source,]
  # print(paste0("Retained ", length(unique(df$locus.id)), " unique loci"))
  # 
  # head(df, n = 3)
  ## /END/ NOT NEEDED? ##
  
  # Per row, tally minor allele variants
  head(df, n = 3)
  print("Per sample, per locus, Tallying variants (homozygous = +2; heterozygous = +1)")
  df$allele.count <- NA
  df[df$Allele.Call=="Homozygous", "allele.count"] <- 2
  df[df$Allele.Call=="Heterozygous", "allele.count"] <- 1
  df[df$Allele.Call=="Absent"|df$Allele.Call=="No Call", "allele.count"] <- 0
  
  # View data
  head(df, n = 10)
  tail(df, n = 10) 
  
  # Build identifier combining the locus.id, locus.hotspot.id and the region.name
  df$locus.id___Region.Name <- paste0(df$locus.id, "___", df$locus.id.rel.hotspot, "___", df$Region.Name)
  head(df)
  
  # Tally minor allele occurrences
  tallies <- df %>% group_by(locus.id___Region.Name) %>%
                summarise(total_views = sum(allele.count))
   
  tallies <- as.data.frame(tallies)
  head(tallies)
  
  # Separate identifier that had combined the locus.id and the region.name
  tallies <- separate(data = tallies, col = "locus.id___Region.Name", into = c("locus.id", "locus.id.rel.hotspot", "Region.Name"), sep = "___", remove = T)
  head(tallies)
  str(tallies)
  nrow(tallies)
  
  # Review the number of possible genotypes at each locus
  head(num_indiv_genod_at_locus.df)
  str(num_indiv_genod_at_locus.df)
  nrow(num_indiv_genod_at_locus.df)

  
  #### 03. Combine alt allele tallies with potential genotypes at locus ####
  nrow(tallies)
  nrow(num_indiv_genod_at_locus.df)
  all_data.df <- merge(x = tallies, y = num_indiv_genod_at_locus.df, by.x = "locus.id.rel.hotspot", by.y = "locus.id", all.x = T)
  nrow(all_data.df)
  head(all_data.df)
  
  # Remove any cases where there was no hotspot value, as these will impact downstream
  all_data.df <- all_data.df[all_data.df$locus.id.rel.hotspot!="NA", ] # note: they are "NA", not NA
  dim(all_data.df)
  length(unique(all_data.df$locus.id.rel.hotspot))
  length(unique(all_data.df$locus.id))
  head(all_data.df)
  
  #### 04. Analyze results ####
  dim(all_data.df)
  
  # Are there any instances of the alt allele being observe more than the total possible? 
  all_data.df[all_data.df$total_views > all_data.df$total.poss.alleles, ]
  nrow(all_data.df[all_data.df$total_views > all_data.df$total.poss.alleles, ])
  
  # Remove these
  all_data.df <- all_data.df[all_data.df$total_views <= all_data.df$total.poss.alleles, ]
  dim(all_data.df)
  
  # Approximate MAF
  all_data.df$MAF <- all_data.df$total_views / all_data.df$total.poss.alleles
  head(all_data.df)
  
  # Ensure MAF is the minor allele
  for(i in 1:nrow(all_data.df)){
    
    #print(i) # Debug only
    
    if(all_data.df$MAF[i] >= 0.5){
      
      all_data.df$MAF[i] <- 1 - all_data.df$MAF[i]
      
    }else if(all_data.df$MAF[i] < 0.5){
      
      all_data.df$MAF[i] <- all_data.df$MAF[i]
      
    }
  }
  
  head(all_data.df, n = 10)
  
  # Export result
  write_delim(x = all_data.df, file = "03_results/tally_all_variants_output.txt", delim = "\t")
  
  # Separate hotspot df and novel df
  hotspot_data.df <- all_data.df[all_data.df$locus.id.rel.hotspot==all_data.df$locus.id, ]
  # Drop NA
  hotspot_data.df$locus.id.rel.hotspot
  dim(hotspot_data.df)
  
  novel_data.df <- all_data.df[all_data.df$locus.id.rel.hotspot!=all_data.df$locus.id, ]
  dim(novel_data.df)
  
  ##### 05. Plotting and comparison #####
  head(hotspot_data.df)
  head(novel_data.df)
  
  novel_data.df$corr.hotspot.MAF <- NA
  
  target <- NULL
  for(i in 1:nrow(novel_data.df)){
    
    target <- novel_data.df$locus.id.rel.hotspot[i]
    
    # What is the corresponding MAF of the hotspot target? 
    novel_data.df$corr.hotspot.MAF[i] <- hotspot_data.df[hotspot_data.df$locus.id.rel.hotspot==target, "MAF"]
    
  }
  
  head(novel_data.df)
  str(novel_data.df)
  
  # Plot
  plot(novel_data.df$MAF, y = novel_data.df$corr.hotspot.MAF)
  
  plot(novel_data.df$corr.hotspot.MAF)
  
  pdf(file = "03_results/test.pdf", width = 100, height = 10)
  boxplot(novel_data.df$MAF ~ novel_data.df$Region.Name.x)
  dev.off()
  
  stripchart(novel_data.df$MAF ~ novel_data.df$Region.Name.x
             , pch = 16
             , cex = 0.5
             )
  
  # Export result
  write_delim(x = novel_data.df, file = "03_results/tally_novel_variants_output.txt", delim = "\t")
  
  #### Key outputs required ####
  # A. Per region, how many alternate variants are there that have higher MAF than the target variant in this dataset? 
  head(novel_data.df)
  
  unique_hotspot_names.vec <- unique(novel_data.df$locus.id.rel.hotspot)
  length(unique_hotspot_names.vec)
  
  target <- NULL; slice <- NULL; slice_data <- NULL; slice_data_all <- NULL
  for(i in 1:length(unique_hotspot_names.vec)){
    
    target <- unique_hotspot_names.vec[i]
    
    slice <- novel_data.df[novel_data.df$locus.id.rel.hotspot==target, ]
    
    num_variants_greater_than_hotspot    <- sum(slice$MAF > slice$corr.hotspot.MAF[1], na.rm = T)
    
    num_variants_greater_than_MAF_cutoff <- sum(slice$MAF > 0.1, na.rm = T)
    
    avg_MAF_in_slice_greater_than_MAF_cutoff <- mean(slice$MAF[slice$MAF > 0.01], na.rm = T)
    
    # Join
    slice_data <- c(target, num_variants_greater_than_hotspot, num_variants_greater_than_MAF_cutoff, avg_MAF_in_slice_greater_than_MAF_cutoff)
    
    slice_data_all <- rbind(slice_data_all, slice_data)
    
  }
  colnames(slice_data_all) <- c("hotspot_ID", "num_var_greater_than_hotspot", "num_var_greater_than_MAF_cutoff", "mean_var_MAF_greater_than_MAF_cutoff")
  slice_data_all <- as.data.frame(slice_data_all)
  
  # Export result
  write_delim(x = slice_data_all, file = "03_results/alternate_novel_variants_summary.txt", delim = "\t")
  
  par(mfrow=c(2,1))
  hist(x = as.numeric(slice_data_all$num_var_greater_than_hotspot), las = 1, main = ""
       , xlab = "Number variants per amplicon with MAF > hotspot MAF")
  hist(x = as.numeric(slice_data_all$num_var_greater_than_MAF_cutoff), las = 1, main = ""
       , xlab = "Number variants per amplicon with MAF > 0.1 (incl. hotspot)")
       
  
  par(mfrow=c(2,1))
  hist(x = as.numeric(slice_data_all$mean_var_MAF_greater_than_MAF_cutoff), las = 1, main = ""
       , xlab = "Average MAF for all variants per amplicon greater than MAF cutoff")
  hist(x = as.numeric(unique(novel_data.df$corr.hotspot.MAF)), las = 1, main = ""
       , xlab = "Per locus MAF for hotspots")  # Note: this uses a workaround that should not be used
  
  # B. What are the chr and position IDs for all the variants? Export as BED
  export_data <- all_data.df$locus.id # note: using all data, so we can get the hotspot as well (still missing a few)
  export_data <- as.data.frame(export_data)
  head(export_data)
  dim(export_data)
  
  export_data  <- separate(data = export_data
                   , col = "export_data"
                   , into = c("chr", "pos")
                   , sep = "__"
                   , remove = T
                    )
          
  # Export result
  write_delim(x = export_data, file = "03_results/all_SNP.txt", delim = "\t")
  
  # Export as BED
  export_data$pos  <- as.numeric(export_data$pos)
  export_data$start <- export_data$pos - 1
  head(export_data)
  export_data.bed <- export_data[,c("chr", "start", "pos")]
  colnames(export_data.bed)[which(colnames(export_data.bed)=="pos")] <- "end"
  head(export_data.bed)
  
  # Export result
  write_delim(x = export_data.bed, file = "03_results/all_SNP.bed", delim = "\t")
  
  
  
}
