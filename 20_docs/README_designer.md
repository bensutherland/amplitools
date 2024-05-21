## Panel designer support ##
This section will use as input a tab-delimited list of marker names and a corresponding VCF file and reference genome to extract a sequence window flanking the target variant and prepare it for submission to a commercial provider for primer panel design.       

#### 01. Obtain information from VCF ####
Use a tab-delimited file with the first column as marker names to extract necessary information from the VCF:     
`./01_scripts/designer/01_collect_info_from_vcf.sh <marker input file>.txt <marker VCF file>.vcf`      

Output:      
`10_designer/vcf_selection.csv`, which will be a comma-delimited file in the format of:     
`NW_018405745.1, 19515, 100055:13:-, C, T`        
This is a text file with fields (1) chr; (2) pos of SNP in ref genome; (3) info about marker; (4) ref allele; (5) alt allele.     

Note: this script will remove any marker from the output that is within 200-bp from the start of the contig, as the full window will not be extractable.       

If you are not starting from a VCF, you can also create your own file with the following format (no spaces), and save this as `10_designer/vcf_selection.csv`:       
<scaffold_name>,<position>,<marker_name>,<ref_allele>,<alt_allele>           
...then move onto the next stage.     

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

