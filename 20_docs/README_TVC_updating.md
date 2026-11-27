# Updating TVC files
In some cases, raw TVC files may need updating of sample names or run names. It is suggested to update this before proceeding with the regular pipeline.           
The following is in-development, please use with caution and report any bugs, thank you!      


#### Requirements: ####
- RStudio

### General Comments ###
This README aims to resolve the following issues, should they arise:       
(1) double underscores (`__`) present in run names; 
(2) changing sample names in the TVC file.   

Item (2) may occur when a sequencing facility puts many other details in the sample name, or if there are mix-ups and corrections to be made to the TVC file after the run has occurred.      

Note: Item (2) requires that Item (1) clean-up function is run.      


### Cleaning up run names ###
Put the TVC file to correct into `02_input_data/tvc_to_update/`. As always, have a back-up of the raw file somewhere else.       
```
# Source amplitools in RStudio to activate functions
clean_run_ids(tvc_file.FN = "") # specify relative path to the TVC file to be fixed

# Output will be in the same folder, with added suffix '_clean_run_id.xls'   

```


### Updating sample IDs in TVC file ###
Once the run name cleanup has been run (see above), do the following:    
```
# Create a fillable sample renaming form including all existing sample names for TVC file(s)
prep_sample_key()
# note: this will look for files ending in '*_clean_run_id.xls' within '02_input_data/tvc_to_update/' 
#  and outputs '02_input_data/tvc_to_update/<tvc_filename>_sample_key.txt'

# Add the target renaming to the sample key and save as '*_clean_run_id_sample_key_annot.txt'   

# Apply the sample ID correction
correct_tvc_ids(sample_key.FN = <path_to_annotated_sample_key>)  
# note: outputs updated TVC file to '02_input_data/<filename>_clean_run_id_updated_ids.xls'
#   also outputs log file of the rename applied '02_input_data/<filename>_clean_run_id_sample_key_annot_applied.txt'   

```

After the above is completed, go to the main README workflow of converting TVC files to genind files.     


[Back to main README](https://github.com/bensutherland/amplitools)    

