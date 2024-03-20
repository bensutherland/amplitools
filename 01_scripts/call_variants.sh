#!/bin/bash
# Call variants from a merged bam

# System variables
INPUT_FOLDER="13_mapped_mhap"
INPUT_FILE="all_merged.bam"
GENOME_FOLDER="00_archive"
GENOME_FILE="cgig_v.1.0_amplicon_ref.fna"
CORES=6
OUTPUT_FOLDER="14_extract_mhap"
OUTPUT_FILE="mpileup_calls.bcf"

# Call variants from merged bam
bcftools mpileup -D -d 115000 \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    -f $GENOME_FOLDER/$GENOME_FILE $INPUT_FOLDER/$INPUT_FILE --threads $CORES | 
    bcftools call -mv --annotate GQ -Ob -o $OUTPUT_FOLDER/$OUTPUT_FILE --threads $CORES

