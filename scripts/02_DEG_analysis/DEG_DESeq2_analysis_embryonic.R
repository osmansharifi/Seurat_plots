library(Seurat)
library(DESeq2)
library(glue)

# By Viktoria Haghani

################################################################################
## Variables

# Paths

#data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_E18_with_labels_proportions.rda"
data_file <- "/share/lasallelab/Osman/scAlign_collab/Male_Cortex_labeled/rett_E18_with_labels_proportions.rda"

#DEG_data_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/M_MUT_and_WT_M_E18_WB/"
DEG_data_dir <- "/share/korflab/home/viki/snRNA-seq-pipeline/DEG_data/DESeq2/M_MUT_and_WT_M_E18_WB/"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_M_E18_WB"

################################################################################
## Data preparation

# Load in data
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'

# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')

## Reorganize Seurat object identities for DEG analysis
# Create only MUT and WT groups
experiment.aggregate@meta.data$new.ident <- plyr::mapvalues(
  x = experiment.aggregate@meta.data$orig.ident, 
  from = c("MUT_M_E18_WB1", "MUT_M_E18_WB2", "WT_M_E18_WB1", "WT_M_E18_WB2"), 
  to = c("MUT_M_E18_WB", "MUT_M_E18_WB", "WT_M_E18_WB", "WT_M_E18_WB")
)

################################################################################
## DESeq2 Analysis

# Add one count to every RNA count so there are no zeroes in data set for DESeq2 log function (pseudocount)
# This is necessary because without pseudocounting, DESeq2 will have an error
experiment.aggregate[["RNA"]]@counts<-as.matrix(experiment.aggregate[["RNA"]]@counts)+1

for (cell_type in cell_types){
  DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_E18_WB", group.by = "new.ident", subset.ident = cell_type, test.use = "DESeq2", slot = "counts")
  cell_DESeq2_DEG <- subset(x = DESeq2_DEG, subset = p_val_adj < 0.05)
  # Write data to CSV so analysis does not need to be rerun when working with data
  write.csv(cell_DESeq2_DEG, file = glue(DEG_data_dir, cell_type, "_", metadata_info, "_DESeq2_DEG.csv"))
}
