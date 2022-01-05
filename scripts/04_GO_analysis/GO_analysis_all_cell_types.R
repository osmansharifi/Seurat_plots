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
time_points <- list("E18_WB", "P30_CORT", "P60_CORT", "P120_CORT")
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo")
go_ontologies <- list("BP", "CC", "MF")

# Variable
metadata <- "M_MUT_and_WT_M"

################################################################################

# Read in GO Data
for (time_point in time_points){
  for (cell_type in cell_types){
    for (ont in go_ontologies){
      assign(glue(cell_type, "_", time_point, "_", ont), read.csv(file = glue(go_data_path, metadata, "_", time_point, "/", cell_type, "_", metadata, "_", time_point, "_", ont, "_gentable.csv")))
    }
  }
}