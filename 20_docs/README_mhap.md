# Microhaplotype Analysis
As with the main analysis, this workflow comes with no guarantees of usefulness. This workflow is currently under development.       

#### Requirements: ####
This workflow should work on either linux or mac OS.      
[bwa](https://github.com/lh3/bwa)       
[samtools](https://samtools.sourceforge.net)      
[bcftools](https://samtools.github.io/bcftools/bcftools.html)       
[bedtools](https://bedtools.readthedocs.io/en/latest/)        
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)      
[multiqc](https://multiqc.info)    
[microTyper](https://github.com/delomast/microTyper)      


### General Comments ###
To identify novel variants, the following approach will be taken:      
(1) create an amplicon-only reference genome to align against (one contig = one mhap);      
(2) align amplicon panel fastq.gz results against the amplicon-only reference;       
(3) call variants from the above, produce a VCF file with all samples;       
(4) inspect VCF file, filter as needed;      
(5) use bam and VCF files to generate microhaplotypes.      


### 01. Generate amplicon-only assembly ###
Downstream application will assume that each contig is a single microhaplotype, and so it is necessary to create a amplicon sequence-only reference.      

Obtain the regions file and reference genome sourced from [amplitargets](https://github.com/bensutherland/amplitargets), and put in `00_archive`.     

Extract the amplicon sequence only and index, as demonstrated for Pacific oyster:     
```
# Decompress the fasta
gunzip 00_archive/GCA_000297895.1_oyster_v9_genomic.fna.gz      

# Create amplicon-only fasta  
bedtools getfasta -fi 00_archive/GCA_000297895.1_oyster_v9_genomic.fna -bed 00_archive/WGAG22008_BJS_OYRv01_Region_A.bed > 00_archive/cgig_v.1.0_amplicon_ref.fna

# index the newly created fasta for downstream alignments
bwa index 00_archive/cgig_v.1.0_amplicon_ref.fna

```       


### 02. Prepare sample data ###
Obtain per-sample demultiplexed fastq data from the AmpliSeq provider, and put this data into the folder `12_input_mhap`.         


Aside: if the per-sample data is in bam format, convert it to fq.gz:     
`01_scripts/bamtofastq.sh`      


### 03. Quality check ###
Use fastqc/ multiqc to evaluate quality:      
```
fastqc 12_input_mhap/*.fq.gz -o 12_input_mhap/fastqc_raw -t 4 
multiqc -o 12_input_mhap/fastqc_raw/ 12_input_mhap/fastqc_raw/    
``` 

### 04. Align samples against reference genome ### 
Update the variable for the reference genome, then run the following script:       
`01_scripts/bwa_mem_align_reads.sh 12`       
Requires that samples suffix is .fastq.    
Note: this will align, sort, index, and generate idxstats. It will add read group IDs to samples, as these are required for downstream analyses.      

idxstats provides a tab-delimited output with each line consisting of a reference sequence name, the sequence length, the number of mapped read segments, and the number of unmapped read-segments.     

Aligned output will be in `13_mapped_mhap`.       


### 05. Call variants from the aligned samples ###
Variants will be called from the aligned samples as follows:      
```
# Prepare a list of all sorted bam files
ls -1 13_mapped_mhap/*.sorted.bam > 13_mapped_mhap/bamlist.txt

# Merge bam files
samtools merge 13_mapped_mhap/all_merged.bam -b 13_mapped_mhap/bamlist.txt --threads 6

# Index the reference genome
samtools faidx 00_archive/cgig_v.1.0_amplicon_ref.fna      

```

Call variants with mpileup by updating any variables and running:       
`./01_scripts/call_variants.sh`     
Output will be in `14_extract_mhap`.         


### 06. Filter the called variants ###
Conduct light filtering:     
`bcftools view -i 'F_missing < 0.1 & TYPE="snp" & QUAL>=20 & FORMAT/DP>10' --min-alleles 2 --max-alleles 2 14_extract_mhap/mpileup_calls.bcf -Ob -o 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1.bcf`

Settings:   
F_missing:      fraction of missing genotypes per locus
TYPE="snp":     keep only SNPs
QUAL:           SNP quality value
DP:             depth (per sample; #TODO: confirm)
--min-alleles:  at least this many alleles observed per locus
--max-alleles:  at most this many alleles observed per locus

How many SNPs remain? 
`bcftools view 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1.bcf | grep -vE '^#' | wc -l`      

Filter on MAF?     
```
bcftools +fill-tags 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1.bcf -Ob -o 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF.bcf  -- -t AF

bcftools view -i 'INFO/AF > 0.01' 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF.bcf -Ob -o 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF_maf0.01.bcf
```


### 07. Call microhaplotypes ###
Here we will use microTyper2.0 to pull the information out of the bam files in the mapped folder, but first we need to create a Position File to provide to microTyper2.0, based on the VCF with all the SNPs in it.         

```
# Create the content of the position file
bcftools view 14_extract_mhap/mpileup_calls_SNP_only_biallelic_q20_dp10_Fmiss_0.1_w_AF_maf0.01.bcf | grep -vE '^#' - | awk '{ print $1 "\t" $2 "\t" "S" "\t" $5 }' - > 14_extract_mhap/position_file_body.txt

# Add a header to the position file
echo -e "Locus \t RefPos \t Type \t ValidAlt" > 14_extract_mhap/position_file.txt && cat 14_extract_mhap/position_file_body.txt >> 14_extract_mhap/position_file.txt

# Clean up
rm 14_extract_mhap/position_file_body.txt
 
```

Before proceeding, make sure that your position file and your reference genome still match in names.    


Now that the position file is made (and ref genome cleaned up if needed), we can extract genotypes as follows:     
```
mtype2 -f 13_mapped_mhap/*.sorted.bam -p 14_extract_mhap/position_file.txt -r 00_archive/cgig_v.1.0_amplicon_ref.fna -d 10 -c .99 -o 14_extract_mhap/genos.txt
```

This will output a file called `14_extract_mhap/genos.txt`.           

Use the R package EFGLmh to convert from long form to wide form, then export to rubias. 

To go from this file to CKMR-sim (#TODO)

