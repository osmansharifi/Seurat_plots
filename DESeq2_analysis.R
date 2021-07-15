library(Seurat)
library(DESeq2)
library(glue)

# Load in data
load('rett_P30_with_labels_proportions.rda')
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'
# Values represent cell numbers for each cell type
before_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident) 

# We want to get rid of the G2M and S phase cells, so subset to keep only G1 cells
experiment.aggregate <- subset(x = experiment.aggregate, subset = cell.cycle == "G1")

# Subset to remove mitochondrial genes
# Set threshold to 0.5%
experiment.aggregate <- subset(x = experiment.aggregate, subset = percent.mito <= "0.5")
# Values represent cell numbers for each cell type
after_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident)

# Reorganize Seurat object identities for DEG analysis
# Create only MUT and WT groups
experiment.aggregate@meta.data$new.ident <- plyr::mapvalues(
  x = experiment.aggregate@meta.data$orig.ident, 
  from = c("MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2"), 
  to = c("MUT_M_P30", "MUT_M_P30", "WT_M_P30", "WT_M_P30")
)

# Add one count to every RNA count so there are no zeroes in data set for DESeq2 log function (pseudocount)
# This is necessary because without pseudocounting, DESeq2 will have an error
experiment.aggregate[["RNA"]]@counts<-as.matrix(experiment.aggregate[["RNA"]]@counts)+1

# Make a list of cell types in the data
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non-neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 
# Run DESeq2 test for every cell type cluster
for(cell_type in cell_types) {
  file_name <- glue('{cell_type}_DESeq2_DEG') 
  file <- glue('{file_name}.csv')
  file_name <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30", group.by = "new.ident", subset.ident = cell_type, test.use = "DESeq2", slot = "counts")
  write.csv(file_name, file = file)
  }