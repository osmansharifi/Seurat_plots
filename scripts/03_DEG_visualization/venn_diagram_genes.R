library(glue)

# By Viktoria Haghani

################################################################################
## Variables

# Paths
EdgeR_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/"
DESeq2_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/"
Limma_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_M_P30_CORT"

################################################################################

for (cell_type in cell_types){
  # Read in data from Limma analysis
  assign(paste0(cell_type, "_Limma_DEG"), read.csv(file = glue(Limma_DEG_dir, cell_type, "_", metadata_info, "_Limma_DEG.csv")))
  # Read in data from EdgeR analysis
  assign(paste0(cell_type, "_EdgeR_DEG"), read.csv(file = glue(EdgeR_DEG_dir, cell_type, "_", metadata_info, "_EdgeR_DEG.csv")))
  # Read in data from DESeq2 analysis
  assign(paste0(cell_type, "_DESeq2_DEG"), read.csv(file = glue(DESeq2_DEG_dir, cell_type, "_", metadata_info, "_DESeq2_DEG.csv")))
}

