# amplitools
Tools for working with amplicon sequencing data.        
In development stage. No guarantees of usefulness.        

Platforms supported:       
- Ion Torrent

Requirements:       
- Linux or Mac operating system
- R (and packages within Rscripts)
- [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)
- [ms_oyster_panel](https://github.com/bensutherland/ms_oyster_panel) for demo analysis
- [CKMRsim](https://github.com/eriqande/CKMRsim)


## Sections ##
[VariantCaller to genepop](#a-variantcaller-to-genepop)          
[Store results](#store-raw-genotypes)                   
[Analyze tech reps](#analyze-tech-reps)              
[Parentage analysis](#parentage-analysis)               
[Genomic coordinates](#genomic-coordinates)              


### A. VariantCaller to genepop ###
This section will allow you to convert the output of the Ion Torrent variantCaller to a genepop file that can be used for analyzing technical replicates or to go into other population genetic or parentage analyses.       

#### 01. Prepare genotype block ####
Put any number of tab-delimited Ion Torrent VariantCaller output (suffix: xls) in the folder `02_input_data`. An example of the data can be found here: [test data](#TOADD)       

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


### B. Store Results ###
This section will provide suggestions as to how to best store amplicon sequence data from variantCaller outputs. (#TODO) 


### C. Analyze technical replicates ###
This section will...      

Open the [simple_pop_stats](#simple_pop_stats) Rscript `01_scripts/simple_pop_stats_start.R`, update the `on_network` variable to FALSE, and then source the script. This will initiate R functions used in this section (#TODO: note: comp tech reps fn remains in dev scripts).      

Use the following command to compare technical replicates from two genepop files:      
```
comp_tech_reps(format_type = "amplitools", max_missing=0.5)      
flags:      
- format_type indicates that the genepops and individual IDs were created in the present repo
- max_missing sets a cutoff where if the missing data per individual is greater than this proportion the sample will be removed and not considered in the technical replicate comparison

This function will also output a genind file to the global environment containing the best of the two technical replicates (and the other samples without technical replicates). Use the following command to save this before clearing space again:      
`save(obj_nr_best, file = "02_input_data/obj_nr_best_2023-05-01.RData")`      


## 03. Population genetic analysis example ##
You can now use the genepop created above, or the genind from (02), to analyze your data in simple_pop_stats or elsewhere. For example, follow the interactive script or adapt similar steps to filter and analyze your genepop as shown in:        
`ms_amplicon_panel/01_scripts/sps_popgen_analysis.R`       


## 04. Estimate parent-offspring and fullsibs using CKMR-Sim




## 03. Other functions ##
Align a FASTA against a reference genome to see where your markers are:       
`./01_scripts/bwa_align.sh`

`samtools view -q 30 05_genome_plot/query_seqs.sam -o 05_genome_plot/query_seqs_above_30_MAPQ.bam`      

`samtools view 05_genome_plot/query_seqs_above_30_MAPQ.bam | awk '{ print $1 }' - | sort -n | uniq -c | sort -nk1 | awk '{ print $1 }' - | uniq -c | sort -nk1 | less`       

`samtools view 05_genome_plot/query_seqs_above_30_MAPQ.bam | sed 's/:://g' | grep 'LR761' - | awk '{ print $1 "\t" $3 "\t" $4 }' - > 05_genome_plot/mnames_to_alt_genome_pos.txt`      

Note: here is where you could technically reduce your list based on markers passing filters (#TODO#).       
Get `retained_loci_keep_targets.txt` from `simple_pop_stats`.    
`awk '{ print "mname_"$1"JH"}' 02_input_data/retained_loci.txt > 02_input_data/retained_loci_keep_targets.txt`
`grep -f 02_input_data/retained_loci_keep_targets.txt 05_genome_plot/mnames_to_alt_genome_pos.txt > 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci.txt`    

`awk '!seen[$1]++' 05_genome_plot/mnames_to_alt_genome_pos.txt > 05_genome_plot/mnames_to_alt_genome_pos_no_dups.txt` 
or
`awk '!seen[$1]++' 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci.txt > 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci_no_dups.txt`


`awk '{ print $2 }' 05_genome_plot/mnames_to_alt_genome_pos_QCd_loci_no_dups.txt | sort -n | uniq -c > 05_genome_plot/markers_no_dups_per_chr.txt`




o integrate:      
Markers included in proton file:     
`awk -F '\t' '{ print $12 }' 02_input_data/R_2022_08_04_09_19_56_user_S5XL-00533-1089-OYR-20220729_7\ \(2\).xls | grep -v 'tvc' | grep -v 'Allele' - | sort | uniq > markers_included_in_proton_file.txt`

