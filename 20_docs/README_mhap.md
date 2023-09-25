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

#### 03. Align to reference genome 
Align reads against your reference genome:     
`01_scripts/bwa_mem_align_reads.sh`       


#### 04. Make list of BAM files ####
`ls -1 13_mapped_mhp/*.sorted.bam > bamlist.txt`        

