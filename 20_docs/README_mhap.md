### Microhaplotype Analysis ###
This analysis will largely use the microhaplot functions and workflows by ngthomas et al. (#cite). 

#### 00. Prepare data
Put raw data in `12_input_data_mhap` folder. If your data is in bam format, do step (1), and if it is already in fq.gz format, go ahead to step 2.     

Index reference genome:    
`bwa index -p <your_prefix> -a is <your_genome>`       

#### 01. Convert from bam to fastq
The input needs to be fastq format, so use the following script to move your data from bam format to fq.gz:      
`01_scripts/bamtofastq.sh`      

#### 02. Quality check
Use fastqc and multiqc to check and visualize quality in your input data:      
```
fastqc 12_input_data_mhp/*.fq.gz -o 12_input_data_mhp/fastqc_raw -t 56
multiqc -o 12_input_data_mhp/fastqc_raw/ 12_input_data_mhp/fastqc_raw/    
``` 

#### 03. Align to reference genome and associated functions #### 
Align reads against your reference genome:     
`01_scripts/bwa_mem_align_reads.sh 48`       
Note: this will align, sort, index, and generate idxstats.     

idxstats provides a tab-delimited output with each line consisting of a reference sequence name, the sequence length, the number of mapped read segments, and the number of unmapped read-segments.     


#### 04. Make list of BAM files ####
`ls -1 13_mapped_mhp/*.sorted.bam > bamlist.txt`        


#### 05. Detect variants using freebayes ####
# Index the reference genome
`samtools faidx <genome>`
(note: must be decompressed)     

# Call variants
Note: requires to install vcflib `https://github.com/vcflib/vcflib`     

Note 2: it was necessary to have proper RG IDs so that the individual samples were scored while SNP finding occurred.     

Note: I installed freebayes using apt-get, but also installed from github to get the freebayes-parallel script. The main freebayes binary did not seem to be included with the github install. 

Once both were installed, run with:     
`freebayes-parallel <(fasta_generate_regions.py ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai 150) 48 -f ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna -L ./bamlist.txt --haplotype-length 0 -kwVa --throw-away-mnps-obs --throw-away-complex-obs > 14_VCF_mhp/sOCP23_noMNP_noComplex_noPriors.vcf`        

This did not seem to work, let's try without the parallel
`freebayes -f ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna -L ./bamlist.txt --haplotype-length 0 -kwVa --throw-away-mnps-obs --throw-away-complex-obs > 14_VCF_mhp/OCP23.vcf`      


