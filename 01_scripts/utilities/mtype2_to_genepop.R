# Convert mtype2 output to genepop format
# Ben Sutherland (VIU)
# 2024-06-25
# note: simplifying names removes preceding directory, and drops .sorted.bam suffix

mtype2_to_genepop <- function(input.FN = "14_extract_mhap/genos.txt"
                              , output.FN = "14_extract_mhap/genos.gen"
                              , simplify_names = TRUE
                              ){
 
  # Set options
  options(tibble.max_extra_cols = 10)
  
  # Load mtype2 output
  mhap.df <- read.delim(file = input.FN, header = T)
  dim(mhap.df)
  head(mhap.df)
  
  # If simplifying names, removes preceding directory, and dropping .sorted.bam suffix
  if(isTRUE(simplify_names)){
    
    print("Simplifying names by removing source directory and '.sorted.bam' suffix.")
    
    # Clean up individual names
    mhap.df$Indiv <- gsub(pattern = ".*/", replacement = "", x = mhap.df$Indiv)
    mhap.df$Indiv <- gsub(pattern = ".sorted.bam", replacement = "", x = mhap.df$Indiv)
    head(mhap.df)
    
  }
  
  # Clean up locus names, as conversion to EFGL format will remove dots, colons, hyphens
  mhap.df$Locus <- gsub(pattern = ":", replacement = "__", x = mhap.df$Locus)
  mhap.df$Locus <- gsub(pattern = "-", replacement = "_", x = mhap.df$Locus)
  
  # Convert from long form to wide form
  mhap.tbl <- as_tibble(mhap.df)
  wideFormat <- mtype2wide(x = mhap.tbl)  
  dim(wideFormat)
  as.data.frame(wideFormat[1:5,1:5])
  wideFormat.df <- as.data.frame(wideFormat)
  rm(wideFormat)
  wideFormat.df[1:5,1:5]
  
  # Convert into EFGLmh format
  mh.dat <- readInData(input = wideFormat.df, pedigreeColumn = "Pop"
                       , nameColumn = "Indiv"
  )
  
  print(paste0("Number of individuals in dataset: ", length(getInds(mh.dat))))
  print(paste0("Number of loci in dataset: ", length(getLoci(mh.dat))))
  
  # Number of alleles per locus
  print("Calculating allele richness")
  allele_rich <- aRich(mh.dat)
  head(allele_rich)
  print("Number of loci with the following number of alleles per locus: ")
  print(table(allele_rich$aRich))
  
  print("Saving out histogram of number alleles per amplicon as 14_extract_mhap/frequency_of_alleles_per_amplicon_hist.pdf.")
  pdf(file = "14_extract_mhap/frequency_of_alleles_per_amplicon_hist.pdf", width = 9, height = 4)
  hist(allele_rich$aRich
       , breaks = 20
       , main = "", las = 1
       , xlab = "Number alleles per amplicon"
  )
  
  #text(x = 25, y = 45, labels = paste0("Total amplicons: ", nrow(allele_rich)))
  dev.off()
  
  # Export to genepop
  print("Exporting to genepop format.")
  exportGenePop(x = mh.dat, filename = output.FN, useIndNames = T)
  
}
