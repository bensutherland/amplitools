#!/bin/bash
# Run from main directory 
# Assumes all fastq have four lines per record
# Assumes fastq files are compressed with gzip

# User set variables
INPUT_FOLDER="12_input_mhap"
SUFFIX=".fastq.gz"

# Move into the directory containing files
cd $INPUT_FOLDER 

# Clear previous runs
rm reads_per_sample.txt 2> /dev/null
rm reads_per_sample_table.txt 2> /dev/null

## Number Reads
# Determine number of reads per file from samples file
for i in $(ls *$SUFFIX) ; 
    do echo $i ;
    echo "$(gunzip -c $i | wc -l)" / 4 | bc ;
    done >> reads_per_sample.txt

# Separate by second line into two columns
sed 'N;s/\n/ /' reads_per_sample.txt > reads_per_sample_table.txt

# Remove intermediate file
rm reads_per_sample.txt

