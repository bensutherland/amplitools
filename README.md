# amplitools
Tools for working with amplicon sequencing data.         
Ben J. G. Sutherland, Ph.D. (Sutherland Bioinformatics).          

**Note**: this software is provided 'as is', without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in action of contract, tort or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.             

The development of this pipeline has been supported by the following organizations: [Support and Funding page](20_docs/funding_support.md).        

#### Platforms supported:       
- Ion Torrent

#### Panels supported:
For a comprehensive list of panels that have been developed or tested through this pipeline, as well as updated and documented associated files, please see the following repository: [amplitargets](https://github.com/bensutherland/amplitargets).      

#### Requirements:       
- Linux or Mac operating system
- devtools       
- R (and packages within Rscripts)
- [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)
- [CKMRsim](https://github.com/eriqande/CKMRsim)
- Demo analysis (optional): [ms_oyster_panel](https://github.com/bensutherland/ms_oyster_panel)       

#### Citation ####
If you find this tool useful, please cite the original article that uses the tool:        
Sutherland et al. 2023. An amplicon panel for high-throughput and low-cost genotyping of Pacific oyster. bioRxiv [2023.08.26.554329](https://www.biorxiv.org/content/10.1101/2023.08.26.554329v1).        

Please also be sure to cite the tools applied within each function.      

## Sections ##
[A. VariantCaller to genepop](#a-variantcaller-to-genepop)          
[B. Analyze tech reps](#b-analyze-technical-replicates)              
[C. Population genetic analyses](#c-population-genetic-analyses)      
[D. Parentage analysis](#d-parentage-analysis)               
[E. Panel designer](#e-panel-designer)      
[F. Database storage](#f-database-storage)                   

## Getting started ##
Clone this repository and change into the main directory.      
```
git clone https://github.com/bensutherland/amplitools.git
cd amplitools

```

Testing out the platform with test data:      
`load_vc(input_folder = "02_input_data", test_only = TRUE)`         
This should load 8024 rows and 14 columns, with 2079 unique markers, from two different samples, and will save the output as input.list for inspection.              


## A. VariantCaller to genepop ##
Convert [VariantCaller](https://www.thermofisher.com/ca/en/home/life-science/sequencing/next-generation-sequencing/ion-torrent-next-generation-sequencing-workflow/ion-torrent-next-generation-sequencing-data-analysis-workflow/ion-torrent-suite-software.html) output to a genepop file for downstream analysis. An example of the input file is provided in [02_input_data/test_data.xls](https://github.com/bensutherland/amplitools/blob/main/02_input_data/test_data.xls)     


#### 00. Prepare inputs and functions
**Important note:** input filename must be alphanumeric characters only connected with hyphens or underscores (i.e., no spaces). Multiple sequential underscores without other alphanumeric characters breaking up the underscores will cause the script to fail (e.g., `this__will_not_work`).      

Copy any number of input files in the folder `02_input_data`.      

Open the Rscript `01_scripts/00_initiator.R` and source the script. This will initiate R functions used in this section.      


#### 01. Load data
**Background**
VariantCaller format interpretation:     
'Absent' means homozygous reference (0101)       
'No Call' means missing data (0000)         
'Heterozygous' means heterozygous (0102)        
'Homozygous' means homozygous variant (0202)        

Please note: all variantCaller input files **must** have been generated using the same hotspot file. The script assumes that all designations of VariantCaller formats are the same for all files.      

##### 01.a. Prepare genotype block ####
In R, use the following function to **select only hotspot markers**, then convert genotype calls to genepop format to output a multilocus genotype matrix (rows: samples; columns: loci):         
`proton_to_genepop(neg_control="BLANK")`          
note: include the exact string for a pattern only within your negative control samples.     

The output will be a tab-delimited text file in `02_input_data/prepped_matrices/`, one per input file. Sample identifiers will be created as <RunName>__<Barcode>__<SampleName>.      

To inspect mean marker coverage per experimental or control sample, see output files:     
`03_results/mean_marker_cov_per_sample_*.txt`          


##### 01.b. Finalize genepop ####
From the terminal, finalize the genepop by running the following script for each genotype block text file:      
`./01_scripts/format_genepop.sh 02_input_data/prepped_matrices/<filename>.txt`      

The output will be a genepop file for each file output as `02_input_data/prepped_genepops/*.gen`       


## B. Analyze technical replicates ##
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

## C. Population genetic analyses ## 
You can now use the genepop created above, or the genind from (02), to analyze your data in `simple_pop_stats` or elsewhere. For example, follow the interactive script or adapt similar steps to filter and analyze your genepop as shown in:        
`ms_amplicon_panel/01_scripts/sps_popgen_analysis.R`       

`simple_pop_stats` is being developed to work well with `amplitools` to facilitate population genetic analyses. As follows are some examples of useful functions, but please see the `simple_pop_stats` repository for a full explanation.       

#### General functions ####
Source `simple_pop_stats`.    

Load the genepop generated by amplitools:      
`load_genepop(datatype = "SNP")`        

Note: this is now possible within amplitools:      
`load_genepops(genepop_folder = "02_input_data/prepped_genepops/", datatype = "SNP", shorten_name = TRUE)`          
...will load all genepops in the genepop folder into a list, each named based on the source genepop name.     

Simplify your individual names from the amplitools format (i.e., `<run>__<barcode>__<indiv>`) to simply the indiv ID:    
`simplify_names(df = obj, format = "amplitools")`     
Note: must use this function for downstream functions to work.    

#### Update sample attributes ####
Create a sample information file to be used for downstream analyses:       
`generate_popmap(df = obj)`     

This will produce a form to complete in `simple_pop_stats/00_archive/my_data_ind-to-pop.txt`      
It will automatically produce a row per sample (indiv), and provide empty fields for pop, an alternate identifier (alt.ID), sex, and whether the sample should be ignored.     
Manually annotate the above file and save it as `my_data_ind-to-pop_annot.txt`.       
Note: do not include any spaces in field names, and do not start populations with numbers.           

Load and population per-sample information based on the completed file above:     
`annotate_from_popmap(df = obj, popmap.FN = "00_archive/my_data_ind-to-pop_annot.txt")`      

See `simple_pop_stats` for other methods, such as filtering, plotting data, analyses, and converting from genepop to rubias format. The below parentage analysis will depend on a rubias file.      


## D. Parentage analysis ##
The parentage analysis is largely dependent on the R package [CKMRsim](https://eriqande.github.io/CKMRsim/) by Eric C. Anderson. If you use the function `ckmr_from_rubias.R`, please cite CKMRsim.        

This section is in development. For now, please see instructions in the associated pipeline: [ms_oyster_panel](https://github.com/bensutherland/ms_oyster_panel).         


#### Plotting ####
It is now possible to graph the relatives based on the output of CKMR-sim, using the following function, with inputs being the paired relatives and associated log likelihood values (.txt format):      
```
graph_relatives(input.FN = "03_results/parent_fs_broodstock_pw_logl_5.txt"   
                            , logl_cutoff = 5      # Additional logl cutoff as needed
                            , drop_string = "G00"  # String constant to remove from sample names
                            , directed = FALSE     # Should the plotted graph be directed? 
                            , plot_width = 5       # Plotting width
                            , plot_height = 5      # Plotting height
                            )
 
Other inputs may include offspring sibship, or parent-offspring relationships.     

``` 


## E. Panel designer ##
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


## F. Database storage ##
This section will provide suggestions as to how to best store amplicon sequence data from variantCaller outputs. (#TODO) 

