## This is a demonstration analysis of VariantCaller input data to accompany the cgig_ms_panel manuscript
# B. Sutherland (2023-08-16)
# Please see the README for cgig_ms_panel for more detailed instructions

#### 01. amplitools input to genepop format ####
# Clear the workspace then source the 00_initiator.R to activate all functions

# Convert all files in 02_input_data from proton to create genepop matrices
proton_to_genepop(hotspot_only=TRUE, neg_control="BLANK")
#  neg_control = name of the negative control wells
#  hotspot_only = hotspot targets only (default and only method as of 2023-08-16)

# Then run `amplitools/01_scripts/format_genepop.sh <filename>` for each file
# You will then have *.gen files prepared in 02_input_data/prepped_matrices/

# Copy the .gen files into simple_pop_stats/02_input_data


#### 02. simple_pop_stats compare tech reps ####
# Clear the workspace and source simple_pop_stats_start.R

## also source comp_tech_reps.R, which is currently in the dev scripts
comp_tech_reps(format_type = "amplitools", max_missing = 0.5)

date <- format(Sys.time(), "%Y-%m-%d")
save(obj_nr_best, file=paste0("02_input_data/obj_nr_best_", date, ".RData"))


#### 03. simple_pop_stats population genetic analysis

# note: the current analysis uses the obj_nr_best object rather than a genepop as input
#   , given the retention of the best tech rep

## the population genetic analysis is all described in full in 
#   'ms_amplicon_panel/01_scripts/sps_popgen_analysis.R' 
#    and will output a rubias formatted file to 'amplitools/03_prepped_data/cgig_all_rubias.txt'


#### 04. amplitools for parentage analysis ####
# Source 00_initiator.R

# Estimate log likelihoods from the data and simulated sibs/ parents, then calculate on your existing data
# Running manually
input.FN = "03_prepped_data/cgig_all_rubias.txt"
cutoff = 5
parent_pop = "VIU_F1"
offspring_pop = "VIU_F2"
#OR
parent_pop = "VIU_F0"
offspring_pop = "VIU_F1"
# then run ckmr_from_rubias.R interactively


# provide logl cutoff
ckmr_from_rubias(input.FN = "03_prepped_data/cgig_all_rubias.txt", parent_pop = "VIU_F1", offspring_pop = "VIU_F2", cutoff = 5)

# 
