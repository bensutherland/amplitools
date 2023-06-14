#!/bin/bash
# Uses collected info from a VCF to create a bed file

# Global variables
INPUT_FOLDER="10_designer"
OUTPUT_FOLDER="10_designer"
INPUT="vcf_selection.csv"
OUTPUT="vcf_selection.bed"

# Prepare bed file
echo "From each position, will subtract 201 bp (start) and add 200 bp as moving from VCF position to bed"

# Keep the chromosome (1), the lower range, upper range, and marker info
awk -F, '{ print $1 "\t" $2-201 "\t" $2+200 "\t" "mname_"$3 }' $INPUT_FOLDER/$INPUT |
    
    # Remove any info after the first colon in the marker name column 
    awk -F: '{ print $1 }' - >  $OUTPUT_FOLDER/$OUTPUT


echo "Completed, your prepared bed file (tab-delim) is in $OUTPUT_FOLDER/$OUTPUT"

