# Microhaplotype Analysis
Note: in development mode only currently (2023-10-31).       

#### Requirements: ####
bwa       
samtools      
bcftools       
bedtools      
fastqc      
multiqc     
microTyper v.2.0      
SNPlift      

**If identifying new variants:**       
vcflib `https://github.com/vcflib/vcflib`     
freebayes (use apt-get)       
also installed from github to get the freebayes-parallel script. The main freebayes binary did not seem to be included with the github install.       


### General Comments ###
There are two approaches that can be taken:       
(1) use the variantCaller SNPs output by TorrentSuite software (ThermoFisher); and         
(2) identify new variants based on alignments against a ref genome.        

Both approaches share some commonalities, and the two approaches are described below, indicated where the two approaches diverge.       

**If need an amplicon-only assembly**      
Obtain the regions file and reference genome sourced from [amplitargets](https://github.com/bensutherland/amplitargets), and put in `00_archive`.     

Demonstrated with Pacific oyster:       
```
getfasta -fi ./00_archive/GCA_000297895.1_oyster_v9_genomic.fna -bed 00_archive/WGAG22008_BJS_OYRv01_Region_A.bed > 00_archive/Cgig_v.1.0_amplicon_ref_2023-10-31.fna

# index the newly created fasta
bwa index 00_archive/Cgig_v.1.0_amplicon_ref_2023-10-31.fna     

# the reference genome for alignment is now the newly created fasta
```       


### 00. Prepare data
The input sequence data for the pipeline will be in either fastq or bam format. Put these files in the folder `12_input_mhap`.         

If the input data is in bam format, proceed to [Step 01](#01-convert-from-bam-to-fq) below to convert to fastq format.        
If the input data is in fastq format, skip Step 01 and proceed to [Step 02](#02-quality-check-sequence-data).       

Note: if you plan to use SNPs that were called by the ThermoFisher variantcaller, per individual put both the compressed VCF files (.vcf.gz) and the accompanying index file (.vcf.gz.tbi) in the folder `12_input_mhap`.      


### 01. Convert from bam to fq
If your input files are in bam format, use the following to convert to fq.gz format:     
`01_scripts/bamtofastq.sh`      

Note: if your input files are already in fastq format, skip to the next step.       


### 02. Quality check sequence data
Now that the data is in fq.gz format, use fastqc/ multiqc to evaluate quality:      
```
fastqc 12_input_mhap/*.fq.gz -o 12_input_mhap/fastqc_raw -t 56
multiqc -o 12_input_mhap/fastqc_raw/ 12_input_mhap/fastqc_raw/    

# note: results will be in the folder 12_input_mhap/fastqc_raw 
``` 


### 03. Alignment 
Index reference genome:    
`bwa index <your_genome>`       

Align reads against your reference genome:     
`01_scripts/bwa_mem_align_reads.sh 48`       
Note: this will align, sort, index, and generate idxstats.     

idxstats provides a tab-delimited output with each line consisting of a reference sequence name, the sequence length, the number of mapped read segments, and the number of unmapped read-segments.     

Important: it was essential to be sure that read group IDs were added to the samples here, as the downstream genotyper will require this. This is done automatically by the above script.      


### 04. Approach A: identify variants using third party tools
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
This path is slower.      


### 05. Approach B: use variants identified by variantCaller
#### 05.a) Combine all individual VCF files
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

#### 05.b) Filter the merged VCF file
Next, use light filtering to remove less reliable variants:        
```
bcftools view -i 'TYPE="snp" & MIN(FMT/DP)>3 & MIN(FMT/GQ)>30' -q 0.01:minor 12_input_mhap/merged.vcf -o 12_input_mhap/merged_filtered.vcf     

# Note: this is only suggested, and the user should determine their own cutoffs.     
# The above keeps a variant if any one individual the following is true:      
#  -min depth > 3 
#  -min geno qual > 30
# ...and if the following is true for the dataset altogether:      
#  -minor allele freq (not non-ref, but true MAF) > 0.01
#  -do not keep multi-nucleotide polymorphisms (MNPs)

```


Optional aside: if you'd like to prove it to yourself, see how many variants are only present in a single VCF:    
`gunzip -c 12_input_mhap/TSVC_variants_IonCode_0103.vcf.gz | bcftools query -f '%CHROM\_%POS\n' - > TSVC_IonCode_0103_variants.txt`         
...and compare to another sample. You will see many variants only present in a single file.      


#### 05.c) Convert merged, filtered VCF file to amplicon-only reference
Next, if you are using a VCF that was oriented to a different reference than that used for alignments, convert your VCF coordinates from the reference A to the coordinates to reference B. This is relevant, for example, if you use VCF oriented towards a full reference genome, but want it oriented towards the amplicon-only genome.      

Use SNPLift to convert VCF positions from the reference genome used for amplicon panel and that used for alignments.      
Clone the repo for SNPLift into the parent folder of the present repo. Change into the SNPLift main directory for the rest of this section.     

Copy the "old" and "new" genomes into `03_genomes`.         
Copy the VCF that is to be converted to orient to the new genome into `04_input_vcf`.      

Ensure that both assemblies are indexed with BWA.    

Update the paths to the old genome, new genome, input VCF, and output VCF in `02_infos/snplift_config.sh`.       
Optional: you may also choose to set `CORRECT_ID=0` so as to prevent the ID column from being recalculated, and therefoer retain your original IDs in the new VCF.     

Run SNPLift:      
`time ./snplift 02_infos/snplift_config.sh`      

The output VCF file should now match with the alignments against the amplicon-only reference genome.      

Copy the new VCF to `amplitools/12_input_mhap/`.       

Exit the SNPLift repo and return to amplitools. 
 

#### 06. Call microhaplotypes ####
Here we will use microTyper2.0 to pull the information out of the bam files in the mapped folder, but first we need to create a Position File to provide to microTyper2.0, based on the VCF with all the SNPs in it.         

```
# Create the content of the position file
grep -vE '^#' 12_input_mhap/merged_filtered_SNPLift_to_amplicons.vcf | awk '{ print $1 "\t" $2 "\t" "S" "\t" $5 }' - > 14_extract_mhap/position_file_body.txt

# Add a header to the position file
sed '1i Locus \t RefPos \t Type \t ValidAlt' 14_extract_mhap/position_file_body.txt > 14_extract_mhap/position_file.txt

# Clean up
rm 14_extract_mhap/position_file_body.txt
 
```

Optional: if there are any characters after the first space in the fasta accessions in the input file, this needs to be cleaned up as well. Otherwise this will cause a mismatch between the ref genome and the position file and bams. Suggested fix, as example:     
`awk '{print $1}' GCA_000297895.1_oyster_v9_genomic.fna > GCA_000297895.1_oyster_v9_genomic_shortname.fna`      


Now that the position file is made (and ref genome cleaned up if needed), we can extract genotypes as follows:     
```
mtype2 -f 13_mapped_mhap/*.bam -p 14_extract_mhap/position_file.txt -r ~/genomes/GCA_000297895.1_oyster_v9_genomic_shortname.fna -d 10 -c .99 -o 14_extract_mhap/genos.txt    
```

This will output a file called `14_extract_mhap/genos.txt`.           


#### Addendum ####


Go to amplitools and obtain the names of all of the alignment files:       
`ls 13_mapped_mhp/*.bam | sed 's/13\_mapped\_mhp\///g' - > 13_mapped_mhp/label.txt`    
"Label file path. This customized label file is a tab-separate file that
contains entries of SAM file name, individual ID, and group label. Required"




