#!/bin/bash
# Use bed file to extract relevant sections of a fasta

# User-set variables
REF_GENOME="GCF_002113885.1_ASM211388v2_genomic.fna"

# Global variables
INPUT_FOLDER="02_input_data"
OUTPUT_FOLDER="10_designer"
BED_FILE="10_designer/vcf_selection.bed"
OUTPUT="vcf_selection.fa"

# Remove any existing output
echo "Removing $OUTPUT to make new file" 
rm $OUTPUT_FOLDER/$OUTPUT

# Extract from fasta, keep the full header information in field 4
bedtools getfasta -fi $INPUT_FOLDER/$REF_GENOME \
    -bed $BED_FILE \
    -name > $OUTPUT_FOLDER/$OUTPUT

echo "Creating a two column tab delimited output to match the selected fasta"

# Make associated csv by separating each fasta accession into a tab delimited file
awk 'BEGIN { RS = ">" ; OFS = "\t" } NR>1 { print $1, $2 }' $OUTPUT_FOLDER/$OUTPUT > $OUTPUT_FOLDER/selected_chr_and_seq.txt


