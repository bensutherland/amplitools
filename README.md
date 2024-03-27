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
[C. Population genetic analyses](#c-population-genetic-analyses)      
[D. Parentage analysis](#d-parentage-analysis)               
[F. Database storage](#f-database-storage)                   

Looking for **microhaplotype workflow**? A separate README is available [here](20_docs/README_mhap.md).       
Looking for **panel design**? A separate README is available [here](20_docs/README_designer.md).       
Original approach for **tech reps**? [here](20_docs/README_tech_reps.md).     

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

## F. Database storage ##
This section provides suggestions to maintain amplicon sequence data from variantCaller outputs for long-term storage.     

We suggest the following approach:      
1. keep all raw data stored separately from amplitools; 
2. copy in a specific raw data runs to a fresh instance of the main repo; 
3. save an output file in standard genetics format (e.g., rubias) as a 'curated dataset'. This will be date-stamped, and then should not be altered. If it needs to be altered, return to the start and re-run.      

Here are several general tips on how to manage this:       

#### Sugg. 1: create a hierarchical folder structure for raw data  ####   
As an example considering Pacific oyster, on the server/ computer, in a safe and permanent location, do the following:      
```
mkdir 00_raw_data
cd 00_raw_data
mkdir cgig
cd cgig
mkdir cgig_v.1.0
```
Once the structure is built, copy raw run files (.xls) into the above folder.      
- within the folder, store a single spreadsheet explaining per chip the **chip name, date run, and projects included**. This will be used to track what is in each raw data file, so specific chips can be then copied into an active amplitools repo.      
- note: it may be best to make the above read-only without sudo priviledges, to be sure that files do not get altered
- also include an md5 alongside each raw run file, with the same name but .md5 appended to the end of the filename.      
- also include a sample interp file alongside each chip (see below).     


#### Sugg. 2: create sample interp file for each chip and store with the raw data chip ####
e.g., `R_<your_chip_name>_sample_interp.txt`        
Include the following columns so this file can also be used as a popmap:      

| indiv | pop | alt.ID | sex | ignore |    
| ----- | --- | ------ | --- | ------ |    
| R_2022_08_04_09_19_56_user_S5XL-00533-1089-OYR-20220729_7_IonCode_0501.fastq | VIU_F2 | F2-03 | NA | NA | 
| R_2022_10_07_13_17_04_user_S5XL-0055-315-Oyster_1_IonCode_0909.fastq | VIU_F1 | BR10  | NA | NA |     
| ... | ... | ... | ... | ... | 

Note: the ignore column will eventually be used so that one can set amplitools to ignore the sample, for example if it is decided that the sample should not be included in building downstream datasets, or perhaps it is the poorer of two technical replicates and should be dropped. This is currently not implemented, but is planned in a future release.        


Importantly, when a new hotspot is used to generate the output, the files will no longer be directly comparable, and so a new folder should be created (e.g., `cgig_v.1.1`) for clarity.     
     
