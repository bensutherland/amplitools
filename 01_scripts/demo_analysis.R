## This is a demonstration analysis of VariantCaller input data to accompany the ms_oyster_panel manuscript
# B. Sutherland (2023-08-16)
# See the README for ms_oyster_panel for more detailed instructions
#  available here: https://github.com/bensutherland/ms_oyster_panel/blob/main/README.md

#### 01. amplitools input to genepop format ####
# Clear the workspace then source the 00_initiator.R to activate all functions

# Convert all files in 02_input_data from proton format to create genepop matrices
proton_to_genepop(neg_control="BLANK")
#  neg_control = name of the negative control wells

# Then run `amplitools/01_scripts/format_genepop.sh <filename>` for each file
# You will then have *.gen files prepared in 02_input_data/prepped_matrices/

# Copy the .gen files into simple_pop_stats/02_input_data


#### 02. simple_pop_stats compare tech reps ####
# Clear the workspace and source simple_pop_stats_start.R

## also source comp_tech_reps.R, which is currently in the dev scripts
#    comp_tech_reps operates on all genepop (*.gen) files in 02_input_data
comp_tech_reps(format_type = "amplitools", max_missing = 0.5)

date <- format(Sys.time(), "%Y-%m-%d")
save(obj_nr_best, file=paste0("02_input_data/obj_nr_best_", date, ".RData"))


#### 03. simple_pop_stats population genetic analysis

# note: the current analysis uses the obj_nr_best object rather than a genepop as input
#   , given the retention of the best tech rep

## the population genetic analysis is all described in full in 
#   'ms_oyster_panel/01_scripts/sps_popgen_analysis.R' 
#    and will output a rubias formatted file to 'amplitools/03_prepped_data/cgig_all_rubias.txt'


#### 04. amplitools for parentage analysis ####
## Run this after running ms_cgig_panel/01_scripts/sps_popgen_analysis.R
# Source 00_initiator.R

# Estimate log likelihoods from the data and simulated sibs/ parents, then calculate on your existing data
# # VIU_F1 vs VIU_F2
# ckmr_from_rubias(input.FN = "03_prepped_data/cgig_all_rubias.txt"
#                  , parent_pop = "VIU_F1"
#                  , offspring_pop = "VIU_F2"
#                  , cutoff = 5
#                  )

# VIU_F1 vs VIU_F2, no monomorph, no multimapper
ckmr_from_rubias(input.FN = "03_prepped_data/cgig_no_monomorphs_no_multimapper.txt"
                 , parent_pop = "VIU_F1"
                 , offspring_pop = "VIU_F2"
                 , cutoff = 5
)




# VIU_F0 vs VIU_F1
# ckmr_from_rubias(input.FN = "03_prepped_data/cgig_all_rubias.txt"
#                  , parent_pop = "VIU_F0"
#                  , offspring_pop = "VIU_F1"
#                  , cutoff = 5
# )

## Optional: 
# In terminal, adjust 03_prepped_data/cgig_all_rubias.txt to update F0 and F1 as a single population to demonstrate
#   the potential issues when including grandparents
## command:
# cat 03_prepped_data/cgig_all_rubias.txt | sed 's/VIU_F0/VIU_parent/g' | sed 's/VIU_F1/VIU_parent/g' > 03_prepped_data/cgig_all_rubias_F0_and_F1_combined_as_VIU_parent.txt

# # VIU_F0 AND VIU_F1 vs VIU_F2
# ckmr_from_rubias(input.FN = "03_prepped_data/cgig_all_rubias_F0_and_F1_combined_as_VIU_parent.txt"
#                  , parent_pop = "VIU_parent"
#                  , offspring_pop = "VIU_F2"
#                  , cutoff = 5
# )

# Generate reports (note: this is now done automatically from the above)
#prep_report(relationship = "PO", offspring_ids = "03_results/offspring_indiv.txt")

# Generate relatedness graphs
graph_relatives(input.FN = "03_results/po_VIU_F1_vs_VIU_F2_pw_logl_5.txt", drop_string = "")
#graph_relatives(input.FN = "03_results/po_VIU_F0_vs_VIU_F1_pw_logl_5.txt", drop_string = "")
graph_relatives(input.FN = "03_results/fs_offsp_VIU_F2_pw_logl_5.txt", drop_string = "")
graph_relatives(input.FN = "03_results/fs_parent_VIU_F1_pw_logl_5.txt", drop_string = "")

