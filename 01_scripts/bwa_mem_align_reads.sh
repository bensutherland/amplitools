#!/bin/bash

# Align reads against reference genome
# Adapted from stacks_workflow by Eric Normandeau 
# https://github.com/enormandeau/stacks_workflow

# Global variables
GENOMEFOLDER="00_archive"
GENOME="cgig_v.1.0_amplicon_ref.fna" # note: this is the prefix given during index
INPUTFOLDER="12_input_mhap"
OUTPUTFOLDER="13_mapped_mhap"
IDXFOLDER="idx_output"
NCPU="$1"

# For all files in the input, align
for file in $(ls -1 "$INPUTFOLDER"/*.fastq)
do
    # Name of uncompressed file
    echo "Aligning file $file"

    name=$(basename "$file")
    ID="@RG\tID:$name\tSM:$name\tPL:IonProton"
    # note: this was updated to be sure to add a unique RG ID for downstream processing

    # Align reads 1 step
    bwa mem -t "$NCPU" -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R "$ID" \
        "$GENOMEFOLDER"/"$GENOME" "$INPUTFOLDER"/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > "$OUTPUTFOLDER"/"${name%.fastq}".bam

    # Samtools sort
    samtools sort --threads "$NCPU" -o "$OUTPUTFOLDER"/"${name%.fastq}".sorted.bam \
        "$OUTPUTFOLDER"/"${name%.fastq}".bam

    # Index output
    samtools index "$OUTPUTFOLDER"/"${name%.fastq}".sorted.bam

    # Generate stats
    samtools idxstats --threads $NCPU "$OUTPUTFOLDER"/"${name%.fastq}".sorted.bam > "$OUTPUTFOLDER"/"$IDXFOLDER"/"${name%.fastq}"_idxstats.txt 

    # Cleanup
    rm "$OUTPUTFOLDER"/"${name%.fastq}".bam

done

