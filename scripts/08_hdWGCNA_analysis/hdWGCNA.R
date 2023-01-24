# This program will perform hdWGCNA analysis on single cell data
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(1234)

# load the snRNA-seq dataset
load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all_male_mouse_cortex.RData")
DimPlot(all_male, group.by='celltype.call', label=TRUE) +
  umap_theme() + NoLegend()

# Prepare Seurat Object for WGCNA
metadata <- all_male@meta.data
timepoint <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[3])
})
all_male@meta.data$Time_Point <- unlist(timepoint)

condition <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[1])
})
all_male@meta.data$Condition <- unlist(condition)

adult_male <- NormalizeData(adult_male, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(adult_male)
adult_male <- ScaleData(adult_male, features = all.genes)
adult_male <- RunPCA(adult_male, features = VariableFeatures(object = adult_male))
adult_male <- RunUMAP(adult_male, dims = 1:10)

# Subset on a value in the object meta data
adult_male <- subset(x = all_male, subset = Time_Point != "E18")

# Set up Seurat object for WGCNA
adult_male <- SetupForWGCNA(
  adult_male,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "postnatal_mouse_cortex" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
adult_male <- MetacellsByGroups(
  seurat_obj = adult_male,
  group.by = c("celltype.call", "Time_Point", "Condition"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype.call' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
adult_male <- NormalizeMetacells(adult_male)

# Set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "INH", # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)