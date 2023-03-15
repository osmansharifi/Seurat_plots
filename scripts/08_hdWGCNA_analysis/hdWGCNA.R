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
load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all.cortex.combined.RData")
setwd("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis")
Idents(all.cortex.combined) <- "celltype.call"
# Set up multithreading
enableWGCNAThreads(nThreads = 8)

# Prepare Seurat Object for WGCNA
metadata <- adult_postnatal@meta.data
timepoint <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[3])
})
adult_postnatal@meta.data$Time_Point <- unlist(timepoint)

genotype <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[1])
})
adult_postnatal@meta.data$genotype <- unlist(genotype)

# Subset on a value in the object meta data
adult_postnatal <- subset(x = all.cortex.combined, subset = Age != "E18")

# Preprocess
adult_postnatal <- NormalizeData(adult_postnatal, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(adult_postnatal)
adult_postnatal <- ScaleData(adult_postnatal)#, features = all.genes)
adult_postnatal <- FindVariableFeatures(adult_postnatal, selection.method = "vst", nfeatures = 2000)
adult_postnatal <- RunPCA(adult_postnatal, features = VariableFeatures(object = adult_postnatal))
adult_postnatal <- RunUMAP(adult_postnatal, dims = 1:20)
DimPlot(adult_postnatal, group.by='celltype.call', label=TRUE) +
  umap_theme() 

# Set up Seurat object for WGCNA
adult_postnatal <- SetupForWGCNA(
  adult_postnatal,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "postnatal_mouse_cortex" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
adult_postnatal <- MetacellsByGroups(
  adult_postnatal,
  group.by = c("celltype.call", "Time_Point", "Sex", "genotype"), # specify the columns in adult_postnatal@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype.call' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
adult_postnatal <- NormalizeMetacells(adult_postnatal)
Idents(adult_postnatal) <- "celltype.call"
# Set up the expression matrix
adult_postnatal <- SetDatExpr(
  adult_postnatal,
  group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Sncg", "Non-neuronal", "Astro", "Oligo", "Lamp5"), # the name of the group of interest in the group.by column
  group.by="celltype.call", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
adult_postnatal <- TestSoftPowers(
  adult_postnatal,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(adult_postnatal)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
ggplot2::ggsave("Softpowerthreshold.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# construct co-expression network:
adult_postnatal <- ConstructNetwork(
  adult_postnatal, soft_power=8,
  setDatExpr=FALSE,
  overwrite_tom = TRUE# name of the topoligical overlap matrix written to disk
)

PlotDendrogram(adult_postnatal, main='hdWGCNA Dendrogram')
ggplot2::ggsave("WGCNA_Dendrogram.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# need to run ScaleData first or else harmony throws an error:
adult_postnatal <- ScaleData(adult_postnatal, features=VariableFeatures(adult_postnatal))

# compute all MEs in the full single-cell dataset
adult_postnatal <- ModuleEigengenes(
  adult_postnatal,
  group.by.vars="orig.ident"
)

# harmonized module eigengenes:
hMEs <- GetMEs(adult_postnatal)

# module eigengenes:
MEs <- GetMEs(adult_postnatal, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
adult_postnatal <- ModuleConnectivity(
  adult_postnatal,
  group.by = 'celltype.call', group_name = c("L2_3_IT", "L4", "L5", "L6", "Sst", "Pvalb", "Vip", "Sncg", "Non-neuronal", "Astro", "Oligo", "Lamp5")
)

# plot genes ranked by kME for each module
p <- PlotKMEs(adult_postnatal, ncol=5)

p
ggplot2::ggsave("kME.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
adult_postnatal <- ModuleExprScore(
  adult_postnatal,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
adult_postnatal <- ModuleExprScore(
  adult_postnatal,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  adult_postnatal,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  adult_postnatal,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE) # depending on Seurat vs UCell for gene scoring


# stitch together with patchwork
wrap_plots(plot_list, ncol=3)

# plot module correlagram
ModuleCorrelogram(adult_postnatal)

# get hMEs from seurat object
MEs <- GetMEs(adult_postnatal, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
adult_postnatal@meta.data <- cbind(adult_postnatal@meta.data, MEs)
# plot with Seurat's DotPlot function
p <- DotPlot(adult_postnatal, features=mods, group.by = 'celltype.call')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

## Compute Correlations
# convert genotype to factor
adult_postnatal$Genotype <- as.factor(adult_postnatal$Genotype)
# convert time point to factor
adult_postnatal$Time_Point <- as.factor(adult_postnatal$Time_Point)
# convert celltype to factor
adult_postnatal$celltype.call <- as.factor(adult_postnatal$celltype.call)
# convert celltype to factor
adult_postnatal$orig.ident <- as.factor(adult_postnatal$orig.ident)

# list of traits to correlate
cur_traits <- c('Genotype', 'Time_Point')

adult_postnatal <- ModuleTraitCorrelation(
  adult_postnatal,
  traits = cur_traits,
  group.by='celltype.call'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(adult_postnatal)

names(mt_cor$cor)
PlotModuleTraitCorrelation(
  adult_postnatal,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 3,
  text_color = 'white',
  high_color = 'purple',
  mid_color = 'black',
  low_color = 'yellow',
  plot_max = 0.1,
  combine=TRUE
)

moduleDendro <- getDendro(MEs, distance = "bicor")
plotDendro(moduleDendro, labelSize = 4, nBreaks = 5)


