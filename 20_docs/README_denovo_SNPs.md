# de novo SNP calling 
As with the main analysis, this workflow comes with no guarantees of usefulness. This workflow is currently under development.       

#### Requirements: ####
This workflow should work on either linux or macOS.      
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)      
[multiqc](https://multiqc.info)    
[bwa](https://github.com/lh3/bwa)       
[samtools](https://samtools.sourceforge.net)      
[bcftools](https://samtools.github.io/bcftools/bcftools.html)       


### General Comments ###
To identify novel variants, the following approach will be taken:      
(1) align amplicon panel fastq.gz results against the amplicon-only reference;       
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

Note: if provided per-sample data is in bam format, convert it to fq.gz:     
`01_scripts/bamtofastq.sh`      

#### Optional: select best replicate ####
To determine the deepest sequenced replicate, determine the number of reads in each fastq file:       
`01_scripts/count_reads.sh`       
...note: adjust variable SUFFIX if your files do not have the suffix '.fastq.gz'.       

Prepare an interpretation file, `00_archive/filename_to_sample_map.txt`, with column headers `filename` (fastq name per file, no path) and `sample_id` (name of sample, identical for replicates).     
To create an interpretation file from scratch, use the following:      
`basename -a 12_input_mhap/*.fastq.gz > 00_archive/filename_to_sample_map.txt`       
 ...then use spreadsheet editor to complete the info.     

In RStudio, source amplitools initiator, then run the following:      
```
select_best_rep_fastq(input_folder = "12_input_mhap", metadata.FN = "00_archive/filename_to_sample_map.txt", counts.FN = "12_input_mhap/reads_per_sample_table.txt")         
```
This will output a histogram of reads per sample, and a reads per sample table with only the best sample retained ('best' is based on number of reads).      


Then in terminal, run the following to move any file that was not set to be retained into the removed files folder:     
```
for file in $(cat 12_input_mhap/remove_files.txt); do mv "$file" 12_input_mhap/removed_files/; done
```

Finally, add any manually-selected files that you would like to keep in the analysis into the `12_input_mhap` folder. If these files are paired-end, put the R2 files into the `removed_files` subfolder, and only use the R1 file to match the amplicon panel output data.    


### 03. Quality check ###
Use fastqc/ multiqc to evaluate quality:      
```
fastqc 12_input_mhap/*.fastq.gz -o 12_input_mhap/fastqc_raw -t 4 
multiqc -o 12_input_mhap/fastqc_raw/ 12_input_mhap/fastqc_raw/    
``` 

### 04. Align samples against reference genome ### 
Update the variable for the reference genome, then run the following script:       
`01_scripts/bwa_mem_align_reads.sh 12`       
Requires that the sample filename suffix is .fastq.gz     
Note: this will align, sort, index, and generate idxstats. It will add read group IDs to samples, as these are required for downstream analyses.      

Note: It is suggested to use a non-compressed reference genome for indexing so the same version can be used for downstream genotyping (cannot be compressed with fastq.gz).    

idxstats provides a tab-delimited output with each line consisting of a reference sequence name, the sequence length, the number of mapped read segments, and the number of unmapped read-segments.     

Aligned output will be in `13_mapped_mhap`.       

To inspect the number of reads and number of alignments per sample:    
`01_scripts/assess_results.sh` will produce summary tables in each folder.     


### 05. Call variants from the aligned samples ###
Variants will be called from the aligned samples as follows:      
```
# Prepare a list of all sorted bam files
ls -1 13_mapped_mhap/*.sorted.bam > 13_mapped_mhap/bamlist.txt

# Merge bam files
samtools merge 13_mapped_mhap/all_merged.bam -b 13_mapped_mhap/bamlist.txt --threads 6

# Index the reference genome
samtools faidx <ref_genome>       

```

Call variants with mpileup by updating any variables and running:       
`./01_scripts/call_variants.sh`     
Output will be in `14_extract_mhap`.         


### 06. Filter the called variants ###
Automated filtering:     
`01_scripts/filter_bcf.sh`     

Filter on MAF?     
```
bcftools +fill-tags 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1.bcf -Ob -o 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF.bcf  -- -t AF

bcftools view -i 'INFO/AF > 0.01' 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF.bcf -Ob -o 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF_maf0.01.bcf
```

