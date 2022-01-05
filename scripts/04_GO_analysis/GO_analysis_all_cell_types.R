library(foreach)
library(glue)
library(ggplot2)

# By Osman Sharifi & Viktoria Haghani

# Note that GO_analysis_enrichment_scores.R should be run before running this script
# This script reads in the GO Terms previously identified

################################################################################
## Variables

#Paths
go_data_path <- "~/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/"
figure_path <- "~/GitHub/snRNA-seq-pipeline/figures/go_analysis/enrichment_scores/"

# Lists
concise_metadatas <- list("M_MUT_and_WT_M_E18_WB", "M_MUT_and_WT_M_P30_CORT", "M_MUT_and_WT_M_P60_CORT", "M_MUT_and_WT_M_P120_CORT")
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo")
ont_types <- list("BP", "CC", "MF")

################################################################################

# Read in GO Data
foreach(met = concise_metadatas, cell_type = cell_types, ont = ont_types) %do% {
  print("placeholder")
  # I want it to look like: L6_P30_CORT
  # Can use assign(glue(variable_name), function)
}