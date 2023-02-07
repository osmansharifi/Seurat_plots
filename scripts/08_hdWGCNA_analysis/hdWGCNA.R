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

# Prepare Seurat Object for WGCNA
metadata <- all_male@meta.data
timepoint <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[3])
})
all_male@meta.data$Time_Point <- unlist(timepoint)

genotype <- lapply(metadata$orig.ident, function(x) {
  split_name <- strsplit(x, "_")[[1]]
  return(split_name[1])
})
all_male@meta.data$genotype <- unlist(genotype)

# Subset on a value in the object meta data
adult_male <- subset(x = all_male, subset = Time_Point != "E18")

# Preprocess
adult_male <- NormalizeData(adult_male, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(adult_male)
adult_male <- ScaleData(adult_male, features = all.genes)
adult_male <- FindVariableFeatures(adult_male, selection.method = "vst", nfeatures = 2000)
adult_male <- RunPCA(adult_male, features = VariableFeatures(object = adult_male))
adult_male <- RunUMAP(adult_male, dims = 1:20)
DimPlot(adult_male, group.by='celltype.call', label=TRUE) +
  umap_theme() 

# Set up Seurat object for WGCNA
adult_male <- SetupForWGCNA(
  adult_male,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "postnatal_mouse_cortex" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
adult_male <- MetacellsByGroups(
  adult_male,
  group.by = c("celltype.call", "Time_Point"), # specify the columns in adult_male@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype.call' # set the Idents of the metacell seurat object
)
## Removing the following groups that did not meet min_cells: Endo#P120#MUT, Endo#P120#WT, Endo#P30#MUT, Endo#P30#WT, Endo#P60#MUT, Endo#P60#WT, Peri#P30#MUT, Peri#P30#WT, Peri#P60#MUT, Peri#P60#WT, Sncg#P30#MUT, Sncg#P30#WT, Sncg#P60#WT ##

# normalize metacell expression matrix:
adult_male <- NormalizeMetacells(adult_male)

# Set up the expression matrix
adult_male <- SetDatExpr(
  adult_male,
  group_name = "P30", # the name of the group of interest in the group.by column
  group.by='Time_Point', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
adult_male <- TestSoftPowers(
  adult_male,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(adult_male)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
adult_male <- ConstructNetwork(
  adult_male, soft_power=5,
  setDatExpr=FALSE,
  overwrite_tom = TRUE# name of the topoligical overlap matrix written to disk
)

PlotDendrogram(adult_male, main='hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
adult_male <- ScaleData(adult_male, features=VariableFeatures(adult_male))

# compute all MEs in the full single-cell dataset
adult_male <- ModuleEigengenes(
  adult_male,
  group.by.vars="orig.ident"
)

# harmonized module eigengenes:
hMEs <- GetMEs(adult_male)

# module eigengenes:
MEs <- GetMEs(adult_male, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
adult_male <- ModuleConnectivity(
  adult_male,
  group.by = 'Time_Point', group_name = "P30"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(adult_male, ncol=5)

p

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
adult_male <- ModuleExprScore(
  adult_male,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
adult_male <- ModuleExprScore(
  adult_male,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  adult_male,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  adult_male,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE) # depending on Seurat vs UCell for gene scoring


# stitch together with patchwork
wrap_plots(plot_list, ncol=3)

# plot module correlagram
ModuleCorrelogram(adult_male)

# get hMEs from seurat object
MEs <- GetMEs(adult_male, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
adult_male@meta.data <- cbind(adult_male@meta.data, MEs)
# plot with Seurat's DotPlot function
p <- DotPlot(adult_male, features=mods, group.by = 'celltype.call')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

## Compute Correlations
# convert genotype to factor
adult_male$Genotype <- as.factor(adult_male$Genotype)
# convert time point to factor
adult_male$Time_Point <- as.factor(adult_male$Time_Point)
# convert celltype to factor
adult_male$celltype.call <- as.factor(adult_male$celltype.call)
# convert celltype to factor
adult_male$orig.ident <- as.factor(adult_male$orig.ident)

# list of traits to correlate
cur_traits <- c('Genotype', 'Time_Point')

adult_male <- ModuleTraitCorrelation(
  adult_male,
  traits = cur_traits,
  group.by='celltype.call'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(adult_male)

names(mt_cor$cor)
PlotModuleTraitCorrelation(
  adult_male,
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

moduleCor <- getCor(MEs, corType = "bicor")
plotHeatmap(moduleCor, rowDendro = moduleDendro, colDendro = moduleDendro)
moduleCorStats <- getMEtraitCor(MEs, colData = MEs, corType = "bicor", robustY = TRUE)

MEtraitCor <- getMEtraitCor(MEs, colData = colData, corType = "bicor")
traitDendro <- getCor(MEs, y = colData, corType = "bicor", robustY = FALSE) %>% getDendro(transpose = TRUE)
plotDendro(traitDendro, labelSize = 3.5, expandY = c(0.65,0.08), file = "Trait_Dendrogram.pdf")

plotMEtraitCor(mt_cor, moduleOrder = moduleDendro$order, traitOrder = traitDendro$order, topOnly = TRUE, label.type = "p", label.size = 4, label.nudge_y = 0, legend.position = c(1.11, 0.795), colColorMargins = c(-1,4.75,0.5,10.1), file = "Top_ME_Trait_Correlation_Heatmap.pdf", width = 8.5, height = 4.25)

plotMEtraitCor(mt_cor, moduleOrder = moduleDendro$order, traitOrder = traitDendro$order)