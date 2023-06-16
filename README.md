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
[VariantCaller to genepop](#variantcaller-to-genepop)          
[Store results](#store-raw-genotypes)                   
[Analyze tech reps](#analyze-tech-reps)              
[Parentage analysis](#parentage-analysis)               
[Genomic coordinates](#genomic-coordinates)              


### VariantCaller to genepop
#### Create genotype block ####
Put any number of tab-delimited Ion Torrent VariantCaller output (*.xls) in `02_input_data`.      
Per file, the following function will convert genotype calls to nucleotide and genepop formats and output a multilocus genotype matrix:        
`proton_to_genepop(hotspot_only=TRUE, neg_control="BLANK")`          
flags:      
- hotspot_only (T/F) will select only the hotspot targets, not novel variants (note: F not implem. yet)
- neg_control is the string indicating negative controls

This will output as `02_input_data/prepped_matrices/*.txt`      

*Notes:*     
Sample identifiers will be created as per `<RunName>__<Barcode>__<SampleName>`.        

VariantCaller format interpretation:     
'Absent' means homozygous reference (0101)       
'No Call' means missing data (0000)         
'Heterozygous' means heterozygous (0102)        
'Homozygous' means homozygous variant (0202)        

Warning: all VC input files must have been generated using the same hotspot file, for one because it assumes that the designated 'reference' and 'variant' alleles for a marker remains constant across all samples.      

#### Finalize genepop ####
For each prepped matrix (see above), run the following bash script to create sample names:      
`./01_scripts/format_genepop.sh <filename>.txt`      
This will output as `02_input_data/prepped_genepops/*.gen`       


## 02. Compare technical replicates ##
Note: if you do not have technical replicates or do not wish to screen them for repeatability, go to Step 3.     

Open `simple_pop_stats_start.R`, clear space, and update `on_network` to FALSE, then source.      

Use the following command to compare technical replicates from two genepop files:      
`comp_tech_reps(format_type = "amplitools", max_missing=0.5)`      
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

