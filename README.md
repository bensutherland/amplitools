# amplitools
Tools for working with amplicon sequencing data.        
In development stage. No guarantees of usefulness.        

#### Platforms supported:       
- Ion Torrent

#### Requirements:       
- Linux or Mac operating system
- R (and packages within Rscripts)
- [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)
- [CKMRsim](https://github.com/eriqande/CKMRsim)
- [ms_oyster_panel](https://github.com/bensutherland/ms_oyster_panel) for demo analysis


## Sections ##
[A. VariantCaller to genepop](#a-variantcaller-to-genepop)          
[B. Store results](#b-result-storage)                   
[C. Analyze tech reps](#c-analyze-technical-replicates)              
[D. Population genetic analyses](#d-population-genetic-analyses)      
[E. Parentage analysis](#e-parentage-analysis)               
[F. Panel designer](#f-panel-designer)      
[G. Genomic coordinates](#g-characterize-genomic-location-of-panel)      

To get started, clone this repo onto your computer and change into the main directory.     

## A. VariantCaller to genepop ##
Convert [VariantCaller](https://www.thermofisher.com/ca/en/home/life-science/sequencing/next-generation-sequencing/ion-torrent-next-generation-sequencing-workflow/ion-torrent-next-generation-sequencing-data-analysis-workflow/ion-torrent-suite-software.html) output to a genepop file for downstream analysis. An example of the input file is provided [here](#toadd)     


#### 00. Data preparation
Ensure that your input filename is only alphanumeric characters connected with hyphens or underscores (no spaces).    

Copy any number of input files in the folder `02_input_data`.      


#### 01. Prepare genotype block ####
Open the Rscript `01_scripts/00_initiator.R` and source the script. This will initiate R functions used in this section.      

In R, use the following function to convert genotype calls to genepop format to output a multilocus genotype matrix (rows: samples; columns: loci):         
```
proton_to_genepop(hotspot_only=TRUE, neg_control="BLANK")          

flags:      
- hotspot_only (T/F; default: T): select only hotspot targets, remove novel variants
- neg_control (string): negative control pattern
Warning: hotspot_only = FALSE is not implemented yet    

The output will be tab-delimited text files in 02_input_data/prepped_matrices/   
Each input variantCaller xls file will be produced as a separate output file.    
Sample identifiers will be created as <RunName>__<Barcode>__<SampleName>.   

VariantCaller format interpretation:     
'Absent' means homozygous reference (0101)       
'No Call' means missing data (0000)         
'Heterozygous' means heterozygous (0102)        
'Homozygous' means homozygous variant (0202)        

```

Please note: all variantCaller input files **must** have been generated using the same hotspot file. The script assumes that all designations of VariantCaller formats are the same for all files.      

#### 02. Finalize genepop ####
Finalize the genepop files by running the following script for each genotype block text file:      
`./01_scripts/format_genepop.sh 02_input_data/prepped_matrices/<filename>.txt`      

The output will be a genepop file for each file output as `02_input_data/prepped_genepops/*.gen`       


## B. Result storage ##
This section will provide suggestions as to how to best store amplicon sequence data from variantCaller outputs. (#TODO) 


## C. Analyze technical replicates ##
This section will...      

Open the [simple_pop_stats](#simple_pop_stats) Rscript `01_scripts/simple_pop_stats_start.R`, update the `on_network` variable to FALSE, and then source the script. This will initiate R functions used in this section (#TODO: note: comp tech reps fn remains in dev scripts).      

Use the following command to compare technical replicates from two genepop files:      
```
comp_tech_reps(format_type = "amplitools", max_missing=0.5)      

# flags:      
# - format_type (string): indicates the source of the genepop and individual IDs (default: amplitools)
# - max_missing (proportion): sets a cutoff for maximum missing genotype proportion to consider a sample for technical replicate comparison

# Note: this will also output a genind to the global environment containing the best individual sample for each unique sample (based on #TODO), where the best is considered the sample with the highest genotyping rate.    
# Retain this object as follows:    
save(obj_nr_best, file = "02_input_data/obj_nr_best_2023-05-01.RData")     

```

## D. Population genetic analyses ## 
You can now use the genepop created above, or the genind from (02), to analyze your data in `simple_pop_stats` or elsewhere. For example, follow the interactive script or adapt similar steps to filter and analyze your genepop as shown in:        
`ms_amplicon_panel/01_scripts/sps_popgen_analysis.R`       


## E. Parentage analysis ##
This section will (#TODO)      


## F. Panel designer ##
This section will use as input a tab-delimited list of marker names and a corresponding VCF file and reference genome to extract a sequence window flanking the target variant and prepare it for submission to a commercial provider for primer panel design.       

#### 01. Obtain information from VCF ####
Use a tab-delimited file with the first column as marker names to extract necessary information from the VCF:     
`./01_scripts/designer/01_collect_info_from_vcf.sh <marker input file>.txt <marker VCF file>.vcf`      

Output:      
`10_designer/vcf_selection.csv`, which will be a comma-delimited file in the format of:     
`NW_018405745.1, 19515, 100055:13:-, C, T`        
This is a text file with fields (1) chr; (2) pos of SNP in ref genome; (3) info about marker; (4) ref allele; (5) alt allele.     

Note: this script will remove any marker from the output that is within 200-bp from the start of the contig, as the full window will not be extractable.       

#### 02. Prepare bed file with selected windows ####
From the output of above, prepare a bed file that has +/- 200 bp positional information from the target SNP:     
`01_scripts/designer/02_prep_bed.sh`         

Output:     
`10_designer/vcf_selection.bed`      

#### 03. Extract window from reference genome ####
Set the variable `REF_GENOME` with the name of the reference genome that has been put in the folder `02_input_data` above, then run the following script to extract the windows around the target variant:      
`01_scripts/designer/03_extract_seq.sh`      

Output:     
`10_designer/vcf_selection.fa`       
`10_designer/selected_chr_and_seq.txt` (a tab-delimited version of the fasta)     

#### 04. Prepare submission form ####
Interactively run the following script to prepare a tab-delimited file that can be used to submit for primer design by a commercial provider:      
`01_scripts/designer/04_make_sub_form.R`        

**Outputs**
- `10_designer/seq_and_minfo_all_data.csv` (full information for data checking)
- `10_designer/seq_and_minfo_for_submission.csv` (submission info only)
The submission csv has fields (1) marker name; (2) chr; (3) ref allele (based on genome); (4) alt allele; (5) strand; (6) marker type; (7) priority level; (8) formatted seq (e.g., ATGC[A/G]ATGC). More details are available in the script.



## G. Characterize genomic location of panel ## 
















