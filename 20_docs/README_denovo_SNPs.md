# de novo SNP calling 
As with the main analysis, this workflow comes with no guarantees of usefulness.        

#### Requirements: ####
- Linux or Mac operating system    
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)      
- [MultiQC](https://multiqc.info)    
- [bwa](https://github.com/lh3/bwa)       
- [SAMtools](https://samtools.sourceforge.net)      
- [bcftools](https://samtools.github.io/bcftools/bcftools.html)       
- [bedtools](https://bedtools.readthedocs.io/en/latest/) for bam to fastq     


### General Comments ###
To identify novel variants, the following approach will be taken:      
(1) align amplicon panel fastq.gz results against the reference genome;       
(2) call variants from the above, produce a VCF file with all samples;       
(3) filter VCF file.      


## Getting started ##
Clone this repository and change into the main directory.      
```
git clone https://github.com/bensutherland/amplitools.git
cd amplitools

```


### 01. Prepare input sequence data ###
Obtain per-sample demultiplexed fastq data from the AmpliSeq provider, and copy it into the folder `12_input_mhap`.         

To save space, compress files if they are not compressed already (required):     
`gzip *.fastq`    

Note: if provided per-sample data is in bam format, convert it to fastq.gz:     
`01_scripts/bamtofastq.sh`      

If there are technical replicates for some samples and you would like some ideas of how to deal with these, see [this page](20_docs/README_select_best_rep.md).      


### 02. Quality check ###
Use fastqc/ multiqc to evaluate quality:      
```
fastqc 12_input_mhap/*.fastq.gz -o 12_input_mhap/fastqc_raw -t 4 
multiqc -o 12_input_mhap/fastqc_raw/ 12_input_mhap/fastqc_raw/    
``` 


### 03. Align samples against reference genome ### 
Update the variable for the reference genome, then run the following script:       
`01_scripts/bwa_mem_align_reads.sh 12`       
Requires that the sample filename suffix is .fastq.gz     
Note: this will align, sort, index, and generate idxstats. It will add read group IDs to samples, as these are required for downstream analyses.      

Note: It is suggested to use a non-compressed reference genome for indexing so the same version can be used for downstream genotyping (cannot be compressed with fastq.gz).    

idxstats provides a tab-delimited output with each line consisting of a reference sequence name, the sequence length, the number of mapped read segments, and the number of unmapped read-segments.     

Aligned output will be in `13_mapped_mhap`.       

To inspect the number of reads and number of alignments per sample:    
`01_scripts/assess_results.sh` will produce summary tables in each folder.     


### 04. Call variants from the aligned samples ###
Variants will be called from the aligned samples as follows:      
```
# Prepare a list of all sorted bam files
ls -1 13_mapped_mhap/*.sorted.bam > 13_mapped_mhap/bamlist.txt

# Merge bam files
samtools merge 13_mapped_mhap/all_merged.bam -b 13_mapped_mhap/bamlist.txt --threads 6

# Index the reference genome
samtools faidx <ref_genome>       

```

Note: if you get an error at the merge stage and have many files, see the [Troubleshooting](README_troubleshooting.md) page.       

Call variants with mpileup by updating any variables and running:       
`./01_scripts/call_variants.sh`     
Output will be in `14_extract_mhap`.         


### 05. Rename samples in BCF file ###
Rename samples from the <run_name><IonCode>.fastq.gz format to a sample name that is usable for the project and will fit with downstream analyses.     
Do this step before filtering, so that you can remove negative control samples (etc) easily, so that these don't affect filters.      
    
```
# Identify the samples in the dataset
bcftools query -l  14_extract_mhap/mpileup_calls.bcf > 14_extract_mhap/samplelist.txt

# Open the sample list text file, add space, then the desired samplename to update
# Save as space-separated, samplelist_rename.txt 
# e.g., 
#oldname newname\n

bcftools reheader --samples 14_extract_mhap/samplelist_rename.txt -o 14_extract_mhap/mpileup_calls_renamed.bcf 14_extract_mhap/mpileup_calls.bcf  

```


### 06. Remove any individuals as needed from the BCF file ###
Remove any samples that are not wanted in the file before filtering the file:     
```
# Identify present samples
bcftools query -l 14_extract_mhap/mpileup_calls_renamed.bcf > 14_extract_mhap/samples_to_retain.txt    

# Delete any rows in the above file manually
# Note: this would also be where you could delete known problematic low GR samples

# Use the file to include only those present after deletions
bcftools view -S 14_extract_mhap/samples_to_retain.txt 14_extract_mhap/mpileup_calls_renamed.bcf -Ob -o 14_extract_mhap/mpileup_calls_renamed_retained.bcf

```


### 07. Filter the called variants ###
Filtering, update the input BCF variable as needed:     
`01_scripts/filter_bcf.sh`     

Note: you will need to change variable for `DATAFOLDER` and `INPUT_BCF`, as well as any of the filtering parameters, as needed.    


Optional: filter by MAF         
```
bcftools +fill-tags 14_extract_mhap/<filtered>.bcf -Ob -o 14_extract_mhap/<filtered>_w_tags.bcf
# Note: the MAF filed will be the minor allele, not just the non-ref allele.     

bcftools view -i 'MAF > 0.05' 14_extract_mhap/*_w_tags.bcf -Ob -o 14_extract_mhap/mpileup_calls_renamed_retained_noindel5_miss0.15_SNP_q99_avgDP10_biallele_minDP10_maxDP10000_minGQ20_miss0.15_w_tags_MAF0.05.bcf
```


### 08. Prepare for downstream analysis ###
Convert the BCF file to a VCF file and then copy to the input folder of the next program (e.g., [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats).          
`bcftools view <input>.bcf -Ov -o <output>.vcf`     


[Back to main README](https://github.com/bensutherland/amplitools)    

