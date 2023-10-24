#!/bin/bash
# Convert bam file to fastq file using bedtools

# System variables
INPUT_FOLDER="12_input_mhap"

# Convert samples from BAM to FASTQ
ls -1 $INPUT_FOLDER/*.bam | 
    sort -u | 
    while read i
    do 
        # Reporting
        echo "Converting $i from bam to fq.gz"

        # Convert from bam to fastq
        bedtools bamtofastq -i $i -fq ${i%.bam}.fq

        # Compress fq
        gzip ${i%.bam}.fq

    done 

