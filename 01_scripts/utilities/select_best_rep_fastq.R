## Select best replicate from multiple fastq files based on numbers of reads
## requires: input fastq files in 12_input_mhap; metadata file connecting filename and samplename
## requires: first run 01_scripts/count_reads.sh to produce 12_input_mhap/reads_per_sample_table.txt 
# B. Sutherland, 2024-06-24

select_best_rep_fastq <- function(input_folder = "12_input_mhap"
                                  , metadata.FN = "00_archive/filename_to_sample_map.txt"
                                  , counts.FN = "12_input_mhap/reads_per_sample_table.txt"
                                 ){
  
  ## Read in data
  # Read in read counts file
  print(paste0("Reading in read counts file from ", counts.FN))
  read_counts.df <- read.delim2(file = counts.FN, sep = " ", header = F)
  colnames(read_counts.df) <- c("filename", "records")
  print(head(read_counts.df))
  print(dim(read_counts.df))
 
  # Read in metadata file (to be used to specify replicates)
  print(paste0("Reading in metadata file from ", metadata.FN))
  metadata.df <- read.delim2(file = metadata.FN, sep = "\t", header = T)
  print(head(metadata.df))
  print(dim(metadata.df))
  
  ## Combine the two files
  # Data checking
  print("Are there any files in your metadata that are not in your input files?")
  print(setdiff(x = metadata.df$filename, y = read_counts.df$filename))
  
  print("Are there any files in your input files that are not in your metadata?")
  print(setdiff(x = read_counts.df$filename, y = metadata.df$filename))
  
  print("Warning: only those files that are in your metadata will be considered for retention")
  
  # Combine
  all_data.df <- merge(x = metadata.df, y = read_counts.df, by = "filename", all.x = T)
  print(head(all_data.df))
  print(dim(all_data.df))
  
  
  ## Sort based on sample_id then number records, descending
  full.df <- all_data.df[with(all_data.df, order(all_data.df$sample_id, all_data.df$records, decreasing = T)), ]
  print(head(full.df))
  
  # Data checking
  print("There are this many rows in your metadata: ")
  print(nrow(full.df))
  print("..and this many unique sample_id: ")
  length(unique(full.df$sample_id))
  
  # Create a df with only one record per unique sample_id
  full.df <- full.df[!duplicated(x = full.df$sample_id), ]
  head(full.df)
  print(paste0("After removing duplicate sample_id, you have ", nrow(full.df), " samples in the dataset."))
  
  ## Creating some summary plots
  print("Creating a summary plot histogram")
  pdf(file = paste0(input_folder, "/reads_per_sample.pdf"), width = 6.5, height = 5)
  hist(all_data.df$records, main = ""
       , breaks = 40, las = 1
       , xlab = "Number reads per sample"
       )
  dev.off()
  
  print(paste0("avg. number reads per sample: ",  round(x = mean(all_data.df$records), digits = 2)))
  print(paste0("median number reads per sample: ",  round(x = median(all_data.df$records), digits = 2)))
  print(paste0("sd number reads per sample: ",  round(x = sd(all_data.df$records), digits = 2)))
  
  ## Write out a non-redundant list of samples that are present in the metadata file
  write_delim(x = full.df, file = paste0(input_folder, "/reads_per_sample_table_no_reps.txt"))
  
  ## Write out a list of samples to not include in downstream analysis based on above (not present in the retain list)
  all_files <- list.files(path = input_folder, pattern = "*.fastq.gz", full.names = T)
  retain_files <- paste0(input_folder, "/", full.df$filename)
  drop_files <- setdiff(x = all_files, y = retain_files)
  drop_files <- as.data.frame(drop_files)
  write.table(x = drop_files, file = paste0(input_folder, "/remove_files.txt"), col.names = F, row.names = F, quote = F)
  
}
