#!/bin/bash
# Run this script from within the amplitools repo for mhap input data analysis
#  assumes single-end data

# Set variables
FQ_SUFFIX=".fastq.gz"
BAM_SUFFIX=".sorted.bam"

# Move into the directory containing input read files
cd 12_input_mhap 

# Clear previous runs
rm reads_per_sample.txt 2> /dev/null
rm reads_per_sample_table.txt 2> /dev/null


## Number Reads
# Determine number of reads per file from samples file
for i in $(ls *$FQ_SUFFIX) ; 
    do echo $i ;
    gunzip -c $i | wc -l | awk 'END { print $1/4}' ;
    done >> reads_per_sample.txt

# Separate by second line into two columns
sed 'N;s/\n/ /' reads_per_sample.txt > reads_per_sample_table.txt

# Remove intermediate file
rm reads_per_sample.txt

cd ..


## Number Mappings
# Move into the directory containing alignment files
cd 13_mapped_mhap 

rm mappings_per_sample.txt 2> /dev/null 
rm mappings_per_sample_table.txt 2> /dev/null

# Determine number of reads per file from samples file
for i in $(ls *$BAM_SUFFIX) ; 
    do echo $i ;
    samtools view $i | wc -l ;
    done >> mappings_per_sample.txt

# Separate by second line into two columns
sed 'N;s/\n/ /' mappings_per_sample.txt > mappings_per_sample_table.txt

# Remove intermediate file
rm mappings_per_sample.txt

cd ..

## Reporting
echo "The assessment is complete." 

