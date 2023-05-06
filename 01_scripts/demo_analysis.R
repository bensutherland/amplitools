# First run 00_initiator.R

# Convert all files in 02_input_data from proton to create genepop matrices
proton_to_genepop(hotspot_only=TRUE, neg_control="BLANK")
# Where neg_control is the name of the negative control wells in the study; and hotspot_only defines if only the hotspot targets would be used

# Then run `amplitools/01_scripts/format_genepop.sh <filename>` for each file
# You will then have *.gen files prepared in 02_input_data/prepped_matrices/

# Next will be in SPS

# Clear workspace and source simple_pop_stats_start.R
## also source comp_tech_reps.R, which is currently in the dev scripts
comp_tech_reps(format_type = "amplitools", max_missing = 0.5)

save(obj_nr_best, file="02_input_data/obj_nr_best_2023-05-01.RData")

# Restart? 
load("02_input_data/obj_nr_best_2023-05-01.RData")

# Next use 
obj_nr_best
# for all downstream analysis

# Complete all described in 'ms_amplicon_panel/01_scripts/sps_popgen_analysis.R'

# Source 00_initiator.R

# Estimate log likelihoods from the data and simulated sibs/ parents, then calculate on your existing data
# provide logl cutoff
ckmr_from_rubias(input.FN = "03_prepped_data/cgig_all_rubias.txt", parent_pop = "VIU_F1", offspring_pop = "VIU_F2", parent_pattern = "BR", cutoff = 5)
