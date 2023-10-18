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
This path was very slow.      


#### 05. Use previously identified variants from flat file ####
Create bed file by using the tally script. 
`01_scripts/dev/tally_to_freq.R`     

This will output `03_results/all_SNP.bed`.       
This path did not work.     


#### 06. Use previously identified variants from VCF ####
To do this, will need to merge all single-sample VCF files provided by the sequencing facility.     

First, create a list of all VCFs to be merged:        
`ls -1 *.vcf.gz > VCF_file_list.txt`      

Combine all VCF files into a single VCF, as follows:      
`bcftools merge --file-list VCF_file_list.txt -o merged.vcf`        

Filtering:        
`vcftools --vcf merged.vcf --max-missing 0.5 --mac 3 --minQ 30 --remove-indels --minDP 3 --recode --recode-INFO-all --out filtered`     

Need to use SNPLift to convert the VCF positions from the reference genome used for amplicon panel and that used for alignments.      
Clone the repo for SNPLift into the parent folder of the present repo. Change into the SNPLift main directory for the rest of this section.     
Copy the chromosome and contig-level assemblies into `03_genomes`.         
Copy the contig-level assembly into `04_input_vcf`.      

Ensure that both assemblies are indexed with BWA.    

Update the following lines in `02_infos/snplift_config.sh`:       

```
export OLD_GENOME="03_genomes/GCA_000297895.1_oyster_v9_genomic.fna"
export NEW_GENOME="03_genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna"

# Output files
export INPUT_FILE="04_input_vcf/filtered.recode.vcf"
export OUTPUT_FILE="contig_to_chr_2023-10-09.vcf"

# Do final corrections to VCF file
export CORRECT_ID=0         # Recompute the ID column from columns 1 and 2 [0, 1].

```
Note: setting the `CORRECT_ID` to 0 above prevents the ID column from being recalculated, so that your original IDs are carried through to the new VCF.       

Run SNPLift:      
`time ./snplift 02_infos/snplift_config.sh`      

The output VCF file provides chromosome locations of the SNP variants.     

Exit the SNPLift repo. 


Go to amplitools and obtain the names of all of the alignment files:       
`ls 13_mapped_mhp/*.bam | sed 's/13\_mapped\_mhp\///g' - > 13_mapped_mhp/label.txt`    
"Label file path. This customized label file is a tab-separate file that
contains entries of SAM file name, individual ID, and group label. Required"




