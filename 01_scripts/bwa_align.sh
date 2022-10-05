#!/bin/bash
# Align the fasta against the reference genome 

# System variables
SAMPLE_FOLDER="02_input_data"
MAPPED_FOLDER="05_genome_plot"

# User-defined variables
THREADS="2"
#GENOME="02_input_data/GCA_000297895.1_oyster_v9_genomic_unwrapped_truncate_header.fa"
GENOME="/Users/wayne/Documents/genomes/Cgig/GCA_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz"

QUERY="vcf_selection_for_Cgig_panel_v.1_panel_design_2022-02-02.fa"

# Align with bwa (requires that ref is indexed) 
bwa mem $GENOME $SAMPLE_FOLDER/$QUERY > $MAPPED_FOLDER/query_seqs.sam

samtools stats $MAPPED_FOLDER/query_seqs.sam > $MAPPED_FOLDER/samstats_query_seqs.txt

