# Expand the region to be used for extracting fasta regions to make a microhap "reference"
#  where each accession should be approximately ~700 bp in size
#  considering limits to the front end, if the contig is too short to expand
# currently implemented for specific purposes, not generalized
# 2023-11-02, Ben J. G. Sutherland

data.bed <- read.delim(file = "00_archive/WGAG22008_BJS_OYRv01_Region_A.bed", skip = 1, header = F)
head(data.bed)

colnames(data.bed) <- c("contig", "start", "end", "region.name", "comment", "gene.id")
head(data.bed)
data.bed$start <- as.numeric(data.bed$start)
data.bed$end <- as.numeric(data.bed$end)

for(i in 1:nrow(data.bed)){
  
  # If the start of the region is 
  if(data.bed$start[i] > 250){
    
    data.bed$start[i] <- data.bed$start[i] - 250
    
  }else { 
    
    data.bed$start[i] <- 0
    
  }
  
  # Add 250 to the end (there is no way to check if that will go over the end of the contig)
  data.bed$end[i] <- data.bed$end[i] + 250
  
}

head(data.bed)
tail(data.bed)

write.table(x = data.bed, file = "00_archive/WGAG22008_BJS_OYRv01_Region_A_expanded.bed"
            , quote = F, col.names = F, row.names = F, sep = "\t")

