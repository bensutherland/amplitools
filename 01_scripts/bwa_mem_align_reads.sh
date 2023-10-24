#!/bin/bash

# Align reads against reference genome
# Adapted from stacks_workflow by Eric Normandeau 
# https://github.com/enormandeau/stacks_workflow

# Global variables
GENOMEFOLDER="/home/greent/genomes"
GENOME="GCA_000297895.1_oyster_v9_genomic.fna" # note: this is the prefix given during index
INPUTFOLDER="12_input_mhap"
OUTPUTFOLDER="13_mapped_mhap"
IDXFOLDER="idx_output"
NCPU="$1"

# For all files in the input, align
for file in $(ls -1 "$INPUTFOLDER"/*.fq.gz)
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
        - > "$OUTPUTFOLDER"/"${name%.fq.gz}".bam

    # Samtools sort
    samtools sort --threads "$NCPU" -o "$OUTPUTFOLDER"/"${name%.fq.gz}".sorted.bam \
        "$OUTPUTFOLDER"/"${name%.fq.gz}".bam

    # Index output
    samtools index "$OUTPUTFOLDER"/"${name%.fq.gz}".sorted.bam

    # Generate stats
    samtools idxstats --threads $NCPU "$OUTPUTFOLDER"/"${name%.fq.gz}".sorted.bam > "$OUTPUTFOLDER"/"$IDXFOLDER"/"${name%.fq.gz}"_idxstats.txt 

    # Cleanup
    rm "$OUTPUTFOLDER"/"${name%.fq.gz}".bam

done
