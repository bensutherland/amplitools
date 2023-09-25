#!/bin/bash

# Align reads against reference genome
# Adapted from stacks_workflow by Eric Normandeau 
# https://github.com/enormandeau/stacks_workflow

# Global variables
GENOMEFOLDER="/home/greent/genomes"
GENOME="GCF_902806645.1_cgigas_uk_roslin_v1_genomic" # note: this is the prefix given during index
INPUTFOLDER="12_input_data_mhp"
OUTPUTFOLDER="13_mapped_mhp"
NCPU="$1"

# For all files in the input, align
for file in $(ls -1 "$INPUTFOLDER"/*.fq.gz)
do
    # Name of uncompressed file
    echo "Aligning file $file"

    name=$(basename "$file")
    ID="@RG\tID:ind\tSM:ind\tPL:IonProton"

    # Align reads 1 step
    bwa mem -t "$NCPU" -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R "$ID" \
        "$GENOMEFOLDER"/"$GENOME" "$INPUTFOLDER"/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > "$OUTPUTFOLDER"/"${name%.fq.gz}".bam


    # Samtools sort
    samtools sort --threads "$NCPU" -o "$OUTPUTFOLDER"/"${name%.fq.gz}".sorted.bam \
        "$OUTPUTFOLDER"/"${name%.fq.gz}".bam

    samtools index "$OUTPUTFOLDER"/"${name%.fq.gz}".sorted.bam

    # Cleanup
    rm "$OUTPUTFOLDER"/"${name%.fq.gz}".bam
done
