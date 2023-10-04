### Microhaplotype Analysis ###
_Note: in development mode only currently (2023-10-04)_       
This analysis will largely use the microhaplot functions and workflows by ngthomas et al. (#cite). 

#### Requirements: ####
bwa       
samtools      


**If identifying new variants**       
vcflib `https://github.com/vcflib/vcflib`     
freebayes using apt-get       
also installed from github to get the freebayes-parallel script. The main freebayes binary did not seem to be included with the github install. 


#### 00. Prepare data
Put raw sequence data in `12_input_data_mhap` folder.        
If the supplied raw sequence data is in bam format, proceed to step (1), otherwise skip to step (2).       

#### 01. Convert bam to fq
The pipeline input is fastq format. Use the following script to convert data from bam to fq.gz:     
`01_scripts/bamtofastq.sh`      


#### 02. Quality check fq data
Use fastqc and multiqc to check and visualize quality in your input data:      
```
fastqc 12_input_data_mhp/*.fq.gz -o 12_input_data_mhp/fastqc_raw -t 56
multiqc -o 12_input_data_mhp/fastqc_raw/ 12_input_data_mhp/fastqc_raw/    
``` 
Results will be in the folder `12_input_data_mhp/fastqc_raw`.       


#### 03. Alignment 
Index reference genome:    
`bwa index -p <your_prefix> -a is <your_genome>`       

Align reads against your reference genome:     
`01_scripts/bwa_mem_align_reads.sh 48`       
Note: this will align, sort, index, and generate idxstats.     

idxstats provides a tab-delimited output with each line consisting of a reference sequence name, the sequence length, the number of mapped read segments, and the number of unmapped read-segments.     

Important: it was essential to be sure that read group IDs were added to the samples here, as the downstream genotyper will require this.      


#### 04. Identify variants
If you already have a VCF that includes all the required sites for extraction of variants within the microhaplotypes, then skip to step (5). If you would like to call variants from your data to produce a new VCF with all observed and filtered variants, complete this step.      

Make a list of all your bam files:        
`ls -1 13_mapped_mhp/*.sorted.bam > bamlist.txt`        

Detect variants using freebayes:      
```
# Index the reference genome (must be decompressed)
samtools faidx <genome> 

# Call variants in parallel
freebayes-parallel <(fasta_generate_regions.py ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai 150) 48 -f ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna -L ./bamlist.txt --haplotype-length 0 -kwVa --throw-away-mnps-obs --throw-away-complex-obs > 14_VCF_mhp/sOCP23_noMNP_noComplex_noPriors.vcf

# OR Call variants in single core mode
freebayes -f ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna -L ./bamlist.txt --haplotype-length 0 -kwVa --throw-away-mnps-obs --throw-away-complex-obs > 14_VCF_mhp/OCP23.vcf 

```











