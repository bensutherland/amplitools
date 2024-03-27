## Analyze technical replicates ##
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
