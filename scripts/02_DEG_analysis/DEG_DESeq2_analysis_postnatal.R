library(Seurat)
library(DESeq2)
library(glue)

# By Viktoria Haghani

################################################################################
## Variables

# Paths

# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_E18_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P30_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P60_with_labels_proportions.rda"
data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P120_with_labels_proportions3.rda"

DEG_data_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_M_P120_CORT"

################################################################################
## Data preparation

# Load in data
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'

## Prepare data
# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')
# We want to get rid of the G2M and S phase cells, so subset to keep only G1 cells
experiment.aggregate <- subset(x = experiment.aggregate, subset = cell.cycle == "G1")
# Set mitochondrial threshold to 0.5%
experiment.aggregate <- subset(x = experiment.aggregate, subset = percent.mito <= "0.5")

## Reorganize Seurat object identities for DEG analysis
# Create only MUT and WT groups
experiment.aggregate@meta.data$new.ident <- plyr::mapvalues(
  x = experiment.aggregate@meta.data$orig.ident, 
  from = c("MUT_M_P120_CORT1", "MUT_M_P120_CORT2", "WT_M_P120_CORT1", "WT_M_P120_CORT2"), 
  to = c("MUT_M_P120_CORT", "MUT_M_P120_CORT", "WT_M_P120_CORT", "WT_M_P120_CORT")
)

################################################################################
## DESeq2 Analysis

# Add one count to every RNA count so there are no zeroes in data set for DESeq2 log function (pseudocount)
# This is necessary because without pseudocounting, DESeq2 will have an error
experiment.aggregate[["RNA"]]@counts<-as.matrix(experiment.aggregate[["RNA"]]@counts)+1

for (cell_type in cell_types){
  DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P120_CORT", group.by = "new.ident", subset.ident = cell_type, test.use = "DESeq2", slot = "counts")
  cell_DESeq2_DEG <- subset(x = DESeq2_DEG, subset = p_val_adj < 0.05)
  # Write data to CSV so analysis does not need to be rerun when working with data
  write.csv(cell_DESeq2_DEG, file = glue(DEG_data_dir, cell_type, "_", metadata_info, "_DESeq2_DEG.csv"))
}
