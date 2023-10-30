### Microhaplotype Analysis ###
_Note: in development mode only currently (2023-10-04)_       

#### Requirements: ####
bwa       
samtools      
bcftools       
fastqc      
multiqc     

_If identifying new variants_       
vcflib `https://github.com/vcflib/vcflib`     
freebayes (use apt-get)       
also installed from github to get the freebayes-parallel script. The main freebayes binary did not seem to be included with the github install.       

#### General Comments ####
There are two approaches that can be taken:       
(1) use the variantCaller SNPs output by TorrentSuite software (ThermoFisher); and         
(2) identify new variants based on alignments against a ref genome.        

Both approaches share some commonalities, and the two approaches are described below, indicated where the two approaches diverge.       


#### 00. Prepare data
The input sequence data for the pipeline will be in either fastq or bam format. Put these files in the folder `12_input_mhap`.         

If the input data is in bam format, proceed to [Step 01](#01-convert-from-bam-to-fq) below to convert to fastq format.        
If the input data is in fastq format, skip Step 01 and proceed to [Step 02](#02-quality-check-sequence-data).       

Note: if you plan to use SNPs that were called by the ThermoFisher variantcaller, per individual put both the compressed VCF files (.vcf.gz) and the accompanying index file (.vcf.gz.tbi) in the folder `12_input_mhap`.      


#### 01. Convert from bam to fq
If your input files are in bam format, use the following to convert to fq.gz format:     
`01_scripts/bamtofastq.sh`      

Note: if your input files are already in fastq format, skip to the next step.       


#### 02. Quality check sequence data
Now that the data is in fq.gz format, use fastqc/ multiqc to evaluate quality:      
```
fastqc 12_input_mhap/*.fq.gz -o 12_input_mhap/fastqc_raw -t 56
multiqc -o 12_input_mhap/fastqc_raw/ 12_input_mhap/fastqc_raw/    

# note: results will be in the folder 12_input_mhap/fastqc_raw 
``` 


#### 03. Alignment 
Index reference genome:    
`bwa index <your_genome>`       

Align reads against your reference genome:     
`01_scripts/bwa_mem_align_reads.sh 48`       
Note: this will align, sort, index, and generate idxstats.     

idxstats provides a tab-delimited output with each line consisting of a reference sequence name, the sequence length, the number of mapped read segments, and the number of unmapped read-segments.     

Important: it was essential to be sure that read group IDs were added to the samples here, as the downstream genotyper will require this. This is done automatically by the above script.      


#### 04. Approach A: identify variants using third party tools
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
This path is very slow.      


#### 05. Approach B: use variants identified by variantCaller
This approach will use SNPs called by the variantCaller (both hotspot targets and novel). The variantCaller produces a VCF per individual.         
Importantly, for novel variants, variantCaller only calls novel sites if the individual's genotype includes a non-ref allele (i.e., heterozygous or homozygous alt).       
Therefore at this stage there will be a lot of missing data, and the filtering below does not filter on missing data yet.          

To use this option, merge all single-sample VCF files:     
```
# Create a list of all VCF files (note: assumes all .vcf.gz are to be merged)
ls -1 12_input_mhap/*.vcf.gz > 12_input_mhap/VCF_file_list.txt

Combine all VCF files into a single VCF, as follows:      
bcftools merge --file-list VCF_file_list.txt -o merged.vcf        
 
```

Next, use light filtering to remove less reliable variants:        
```
bcftools view --types snps -i 'MIN(FMT/DP)>3 & MIN(FMT/GQ)>30' -q 0.01:minor 12_input_mhap/merged.vcf -o 12_input_mhap/merged_filtered.vcf

# Note: this is only suggested, and the user should determine their own cutoffs.     
# The above keeps a variant if any one individual the following is true:      
#  -min depth > 3 
#  -min geno qual > 30
# ...and if the following is true for the dataset altogether:      
#  -minor allele freq (not non-ref, but true MAF) > 0.01

```


Optional aside: if you'd like to prove it to yourself, see how many variants are only present in a single VCF:    
`gunzip -c 12_input_mhap/TSVC_variants_IonCode_0103.vcf.gz | bcftools query -f '%CHROM\_%POS\n' - > TSVC_IonCode_0103_variants.txt`         
...and compare to another sample. You will see many variants only present in a single file.      


#### 06. Call microhaplotypes ####





#### Addendum ####
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




