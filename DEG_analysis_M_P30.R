library(Seurat)
library(DESeq2)
library(limma)
library(MAST)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(BioVenn)
library(biomaRt)


################################################################################
# Visualization and data preparation
# By Osman Sharifi & Viktoria Haghani

## Load in data
load('rett_P30_with_labels_proportions.rda')
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'
# Values represent cell numbers for each cell type
before_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident) 

## Subset cells in G1 and visualize UMAP
# Visualize UMAP
pcs.use <- 10
experiment.aggregate <- RunUMAP(experiment.aggregate, dims = 1:pcs.use)
# This will show us total cells including those in G2M and S phase
DimPlot(experiment.aggregate, reduction = "umap", group.by = "cell.cycle") +
  ggtitle("Cell Type Grouping Including G2M and S Phase Cells")
# We want to get rid of the G2M and S phase cells, so subset to keep only G1 cells
experiment.aggregate <- subset(x = experiment.aggregate, subset = cell.cycle == "G1")
# Generating a UMAP plot to validate that G2M and S phase cells were removed
DimPlot(experiment.aggregate, reduction = "umap", group.by = "cell.cycle") +
  ggtitle("Cell Type Grouping for Subsetted Data (G1 Only)")
# Validate removal of G1 using phase
DimPlot(experiment.aggregate, reduction = "umap", group.by = "Phase") +
  ggtitle("Cell Type Grouping for Subsetted Data from Phase")
# Generate UMAP plot with cell types
DimPlot(experiment.aggregate, reduction = "umap", group.by = "celltype.call", label = TRUE) +
  ggtitle("Cell Types after G1 Subsetting")

## Visualize PCA for each experimental and control group (same as UMAP)
DimPlot(experiment.aggregate, dims = c(1,2), group.by = "orig.ident") +
  ggtitle("Cell Types by Experimental and Control Groups after G1 Subsetting")
# Visualize PCA for cell types
DimPlot(experiment.aggregate, dims = c(1,2), group.by = "celltype.call", label = TRUE) +
  ggtitle("Cell Types After G1 Subsetting")

## Subset to remove mitochondrial genes
# See percent mitochondrial genes before subsetting to threshold
VlnPlot(experiment.aggregate, features = "percent.mito") +
  ggtitle("% mt gene expression before subsetting")
# Set threshold to 0.5%
experiment.aggregate <- subset(x = experiment.aggregate, subset = percent.mito <= "0.5")
# Validate the mitochondrial genes are removed
VlnPlot(experiment.aggregate, features = "percent.mito") +
  ggtitle("% mt gene expression after setting 0.5%  threshold")
# Values represent cell numbers for each cell type
after_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident)

## See changes in cell number after subsetting
# ?all table to compare cells before G2M, S, and mt exclusion
before_subset_cell_counts
# Generate a data frame from the table
before_subset_cell_counts_df <- data.frame(before_subset_cell_counts)
# Order values so bars appear in descending value order
before_subset_cell_counts_df$Var1 <- reorder(before_subset_cell_counts_df$Var1,-before_subset_cell_counts_df$Freq)
sample_name <- before_subset_cell_counts_df$Var2
# Create grouped bar plot
ggplot(before_subset_cell_counts_df, aes(fill=sample_name, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Cell Counts Per Condition Before Subsetting") +
  xlab("Cell Types") +
  ylab("Number of Cells") +
  ylim(0, 3000) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# ?all table to compare cells left after  G2M, S, and mt exclusion
after_subset_cell_counts 
# Generate a data frame from the table
after_subset_cell_counts_df <- data.frame(after_subset_cell_counts)
# Order values so bars appear in descending value order
after_subset_cell_counts_df$Var1 <- reorder(after_subset_cell_counts_df$Var1,-after_subset_cell_counts_df$Freq)
sample_name <- after_subset_cell_counts_df$Var2
# Create grouped bar plot
ggplot(after_subset_cell_counts_df, aes(fill=sample_name, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Cell Counts Per Condition After Subsetting") +
  xlab("Cell Types") +
  ylab("Number of Cells") +
  ylim(0, 3000) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

## Visualize clusters after G2M, S, and mt exclusion
# Generate UMAP plot with cell types after G1 and mitochondrial subsetting
DimPlot(experiment.aggregate, reduction = "umap", group.by = "celltype.call", label = FALSE) +
  ggtitle("Cell Types After G1 & mt Subsetting")
  ggplot2::ggsave("Cell Types After G1 & mt Subsetting.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# Visualize PCA for cell types after G1 and mitochondrial subsetting
DimPlot(experiment.aggregate, dims = c(1,2), group.by = "celltype.call", label = TRUE) +
  ggtitle("Cell Types After G1 & mt Subsetting")
# Visualize UMAP plot grouped by experimental adn control groups to see how cell types match each group
DimPlot(experiment.aggregate, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Cell Types per Group After G1 & mt Subsetting")

## Reorganize Seurat object identities for DEG analysis
# Create only MUT and WT groups
experiment.aggregate@meta.data$new.ident <- plyr::mapvalues(
  x = experiment.aggregate@meta.data$orig.ident, 
  from = c("MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2"), 
  to = c("MUT_M_P30_CORT", "MUT_M_P30_CORT", "WT_M_P30_CORT", "WT_M_P30_CORT")
)

## See counts for orig.ident and validate that they're combined correctly for new.idents
old_ident_counts <- table(experiment.aggregate@meta.data$orig.ident)
old_ident_counts
new_ident_counts <- table(experiment.aggregate@meta.data$new.ident)
new_ident_counts
# See grouping before
DimPlot(experiment.aggregate, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Cell Types per Group Before Combining MUT and WT")
# See grouping after to validate proper separation and data structure
DimPlot(experiment.aggregate, reduction = "umap", group.by = "new.ident") +
  ggtitle("Cell Types per Group After Combining MUT and WT")

################################################################################
# DEG analysis using tests from Seurat (Wilcoxon and MAST)
# By Viktoria Haghani

###### For loop to be made ###### 
# Make a list of all DE tests being used through Seurat
# Can also use "bimod", "roc", "t", "poisson", "negbinom", and "LR"
tests <- list("wilcox", "MAST")
# Make a list of cell types in the data
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 
# Make a list to store names of the data generated
DEG_data <- list()
# Run every test for every cell type cluster
for(test in tests) {
  for(cell_type in cell_types)
  {
    print(test)
    print(cell_type)
  }}

L2_3_IT_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L2_3_IT", test.use = "wilcox")
write.csv(L2_3_IT_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_wilcox_DEG_all_genes.csv")
L2_3_IT_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_wilcox_DEG_all_genes.csv")
L2_3_IT_wilcox_DEG <- subset(x = L2_3_IT_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L2_3_IT_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_wilcox_DEG_only_stat_sig.csv")
L2_3_IT_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L2_3_IT_wilcox_DEG_stat_sig

L6_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L6", test.use = "wilcox")
write.csv(L6_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_wilcox_DEG_all_genes.csv")
L6_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_wilcox_DEG_all_genes.csv")
L6_wilcox_DEG <- subset(x = L6_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L6_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_wilcox_DEG_only_stat_sig.csv")
L6_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L6_wilcox_DEG_stat_sig

Sst_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sst", test.use = "wilcox")
write.csv(Sst_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_wilcox_DEG_all_genes.csv")
Sst_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_wilcox_DEG_all_genes.csv")
Sst_wilcox_DEG <- subset(x = Sst_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Sst_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_wilcox_DEG_only_stat_sig.csv")
Sst_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Sst_wilcox_DEG_stat_sig

L5_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L5", test.use = "wilcox")
write.csv(L5_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_wilcox_DEG_all_genes.csv")
L5_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_wilcox_DEG_all_genes.csv")
L5_wilcox_DEG <- subset(x = L5_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L5_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_wilcox_DEG_only_stat_sig.csv")
L5_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L5_wilcox_DEG_stat_sig

L4_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L4", test.use = "wilcox")
write.csv(L4_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_wilcox_DEG_all_genes.csv")
L4_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_wilcox_DEG_all_genes.csv")
L4_wilcox_DEG <- subset(x = L4_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L4_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_wilcox_DEG_only_stat_sig.csv")
L4_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L4_wilcox_DEG_stat_sig

Pvalb_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Pvalb", test.use = "wilcox")
write.csv(Pvalb_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_wilcox_DEG_all_genes.csv")
Pvalb_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_wilcox_DEG_all_genes.csv")
Pvalb_wilcox_DEG <- subset(x = Pvalb_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Pvalb_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_wilcox_DEG_only_stat_sig.csv")
Pvalb_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Pvalb_wilcox_DEG_stat_sig

Sncg_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sncg", test.use = "wilcox")
write.csv(Sncg_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_wilcox_DEG_all_genes.csv")
Sncg_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_wilcox_DEG_all_genes.csv")
Sncg_wilcox_DEG <- subset(x = Sncg_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Sncg_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_wilcox_DEG_only_stat_sig.csv")
Sncg_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Sncg_wilcox_DEG_stat_sig

Non_neuronal_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Non_neuronal", test.use = "wilcox")
write.csv(Non_neuronal_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_wilcox_DEG_all_genes.csv")
Non_neuronal_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_wilcox_DEG_all_genes.csv")
Non_neuronal_wilcox_DEG <- subset(x = Non_neuronal_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Non_neuronal_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_wilcox_DEG_only_stat_sig.csv")
Non_neuronal_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Non_neuronal_wilcox_DEG_stat_sig

Oligo_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Oligo", test.use = "wilcox")
write.csv(Oligo_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_wilcox_DEG_all_genes.csv")
Oligo_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_wilcox_DEG_all_genes.csv")
Oligo_wilcox_DEG <- subset(x = Oligo_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Oligo_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_wilcox_DEG_only_stat_sig.csv")
Oligo_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Oligo_wilcox_DEG_stat_sig

Vip_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Vip", test.use = "wilcox")
write.csv(Vip_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_wilcox_DEG_all_genes.csv")
Vip_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_wilcox_DEG_all_genes.csv")
Vip_wilcox_DEG <- subset(x = Vip_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Vip_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_wilcox_DEG_only_stat_sig.csv")
Vip_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Vip_wilcox_DEG_stat_sig

Lamp5_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Lamp5", test.use = "wilcox")
write.csv(Lamp5_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_wilcox_DEG_all_genes.csv")
Lamp5_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_wilcox_DEG_all_genes.csv")
Lamp5_wilcox_DEG <- subset(x = Lamp5_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Lamp5_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_wilcox_DEG_only_stat_sig.csv")
Lamp5_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Lamp5_wilcox_DEG_stat_sig

Astro_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Astro", test.use = "wilcox")
write.csv(Astro_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_wilcox_DEG_all_genes.csv")
Astro_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_wilcox_DEG_all_genes.csv")
Astro_wilcox_DEG <- subset(x = Astro_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Astro_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_wilcox_DEG_only_stat_sig.csv")
Astro_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Astro_wilcox_DEG_stat_sig

Peri_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Peri", test.use = "wilcox")
write.csv(Peri_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_wilcox_DEG_all_genes.csv")
Peri_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_wilcox_DEG_all_genes.csv")
Peri_wilcox_DEG <- subset(x = Peri_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Peri_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_wilcox_DEG_only_stat_sig.csv")
Peri_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Peri_wilcox_DEG_stat_sig

Endo_wilcox_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Endo", test.use = "wilcox")
write.csv(Endo_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_wilcox_DEG_all_genes.csv")
Endo_wilcox_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_wilcox_DEG_all_genes.csv")
Endo_wilcox_DEG <- subset(x = Endo_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Endo_wilcox_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_wilcox_DEG_only_stat_sig.csv")
Endo_wilcox_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_wilcox_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Endo_wilcox_DEG_stat_sig

L2_3_IT_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L2_3_IT", test.use = "MAST")
write.csv(L2_3_IT_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_MAST_DEG_all_genes.csv")
L2_3_IT_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_MAST_DEG_all_genes.csv")
L2_3_IT_MAST_DEG <- subset(x = L2_3_IT_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L2_3_IT_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_MAST_DEG_only_stat_sig.csv")
L2_3_IT_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L2_3_IT_MAST_DEG_stat_sig

L6_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L6", test.use = "MAST")
write.csv(L6_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_MAST_DEG_all_genes.csv")
L6_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_MAST_DEG_all_genes.csv")
L6_MAST_DEG <- subset(x = L6_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L6_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_MAST_DEG_only_stat_sig.csv")
L6_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L6_MAST_DEG_stat_sig

Sst_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sst", test.use = "MAST")
write.csv(Sst_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_MAST_DEG_all_genes.csv")
Sst_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_MAST_DEG_all_genes.csv")
Sst_MAST_DEG <- subset(x = Sst_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Sst_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_MAST_DEG_only_stat_sig.csv")
Sst_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Sst_MAST_DEG_stat_sig

L5_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L5", test.use = "MAST")
write.csv(L5_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_MAST_DEG_all_genes.csv")
L5_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_MAST_DEG_all_genes.csv")
L5_MAST_DEG <- subset(x = L5_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L5_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_MAST_DEG_only_stat_sig.csv")
L5_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L5_MAST_DEG_stat_sig

L4_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L4", test.use = "MAST")
write.csv(L4_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_MAST_DEG_all_genes.csv")
L4_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_MAST_DEG_all_genes.csv")
L4_MAST_DEG <- subset(x = L4_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L4_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_MAST_DEG_only_stat_sig.csv")
L4_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- L4_MAST_DEG_stat_sig

Pvalb_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Pvalb", test.use = "MAST")
write.csv(Pvalb_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_MAST_DEG_all_genes.csv")
Pvalb_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_MAST_DEG_all_genes.csv")
Pvalb_MAST_DEG <- subset(x = Pvalb_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Pvalb_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_MAST_DEG_only_stat_sig.csv")
Pvalb_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Pvalb_MAST_DEG_stat_sig

Sncg_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sncg", test.use = "MAST")
write.csv(Sncg_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_MAST_DEG_all_genes.csv")
Sncg_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_MAST_DEG_all_genes.csv")
Sncg_MAST_DEG <- subset(x = Sncg_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Sncg_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_MAST_DEG_only_stat_sig.csv")
Sncg_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Sncg_MAST_DEG_stat_sig

Non_neuronal_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Non_neuronal", test.use = "MAST")
write.csv(Non_neuronal_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_MAST_DEG_all_genes.csv")
Non_neuronal_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_MAST_DEG_all_genes.csv")
Non_neuronal_MAST_DEG <- subset(x = Non_neuronal_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Non_neuronal_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_MAST_DEG_only_stat_sig.csv")
Non_neuronal_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Non_neuronal_MAST_DEG_stat_sig

Oligo_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Oligo", test.use = "MAST")
write.csv(Oligo_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_MAST_DEG_all_genes.csv")
Oligo_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_MAST_DEG_all_genes.csv")
Oligo_MAST_DEG <- subset(x = Oligo_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Oligo_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_MAST_DEG_only_stat_sig.csv")
Oligo_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Oligo_MAST_DEG_stat_sig

Vip_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Vip", test.use = "MAST")
write.csv(Vip_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_MAST_DEG_all_genes.csv")
Vip_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_MAST_DEG_all_genes.csv")
Vip_MAST_DEG <- subset(x = Vip_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Vip_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_MAST_DEG_only_stat_sig.csv")
Vip_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Vip_MAST_DEG_stat_sig

Lamp5_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Lamp5", test.use = "MAST")
write.csv(Lamp5_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_MAST_DEG_all_genes.csv")
Lamp5_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_MAST_DEG_all_genes.csv")
Lamp5_MAST_DEG <- subset(x = Lamp5_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Lamp5_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_MAST_DEG_only_stat_sig.csv")
Lamp5_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Lamp5_MAST_DEG_stat_sig

Astro_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Astro", test.use = "MAST")
write.csv(Astro_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_MAST_DEG_all_genes.csv")
Astro_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_MAST_DEG_all_genes.csv")
Astro_MAST_DEG <- subset(x = Astro_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Astro_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_MAST_DEG_only_stat_sig.csv")
Astro_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Astro_MAST_DEG_stat_sig

Peri_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Peri", test.use = "MAST")
write.csv(Peri_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_MAST_DEG_all_genes.csv")
Peri_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_MAST_DEG_all_genes.csv")
Peri_MAST_DEG <- subset(x = Peri_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Peri_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_MAST_DEG_only_stat_sig.csv")
Peri_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Peri_MAST_DEG_stat_sig

Endo_MAST_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Endo", test.use = "MAST")
write.csv(Endo_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_MAST_DEG_all_genes.csv")
Endo_MAST_DEG <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_MAST_DEG_all_genes.csv")
Endo_MAST_DEG <- subset(x = Endo_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Endo_MAST_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_MAST_DEG_only_stat_sig.csv")
Endo_MAST_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_MAST_DEG_only_stat_sig.csv")
DEG_data[[length(DEG_data) + 1]] <- Endo_MAST_DEG_stat_sig

# Make a heatmap with rows and columns separated by a predefined label using ComplexHeatmap
plotDEG = function(data,       # Data matrix
                   row.labels, # Row annotation for data, assumed to be a vector of characters
                   col.labels, # Column annotation for data, assumed to be a vector of characters
                   filename="~/DEG_sample_code.png"){
  column_anno = HeatmapAnnotation(celltype=as.character(col.labels))
  ## Plot proportion heatmap
  DEG.heatmap = Heatmap(data,
                        top_annotation=column_anno,
                        col = c('purple', 'yellow'),
                        cluster_rows = T,
                        cluster_columns = T,
                        show_row_names = T,
                        row_split=factor(row.labels),
                        column_split=factor(col.labels),
                        show_heatmap_legend=T,
                        border=T)
  ## Plot
  png(filename, width=12, height=12, units='in', res=300)
  draw(DEG.heatmap)
  dev.off()
}

## Make heatmap for Wilcoxon analysis
# Integrate individual data frames to create a DEG heatmap
data_frame_list_wilcox_DEG = list(L2_3_IT_wilcox_DEG, 
                                  L6_wilcox_DEG,
                                  Sst_wilcox_DEG,
                                  L5_wilcox_DEG,
                                  L4_wilcox_DEG,
                                  Pvalb_wilcox_DEG,
                                  Sncg_wilcox_DEG,
                                  Non_neuronal_wilcox_DEG,
                                  Oligo_wilcox_DEG,
                                  Vip_wilcox_DEG,
                                  Lamp5_wilcox_DEG,
                                  Astro_wilcox_DEG,
                                  Peri_wilcox_DEG,
                                  Endo_wilcox_DEG)





for (i in 1:14){
  print(data_frame_list_wilcox_DEG[[i]][["X"]])
}

lapply(data_frame_list_wilcox_DEG, rownames)
gene_set_wilcox

gene_set_wilcox = unique(unlist(lapply(data_frame_list_wilcox_DEG, rownames)))
ct_set_wilcox = unique(colnames(Reduce(cbind, data_frame_list_wilcox_DEG)))
DEG_matrix_wilcox = matrix(0, nrow=length(gene_set_wilcox), ncol=length(ct_set_wilcox), dimnames-list(gene_set_wilcox, ct_set_wilcox))
for(df in data_frame_list_wilcox_DEG){
  DEG_matrix_wilcox[rownames(df), colnames(df)] = df[,1]
}
DEG_df_wilcox = as.data.frame(DEG_matrix_wilcox)
plotDEG = function(DEG_df_wilcox,
                   gene_set_wilcox,
                   ct_set_wilcox,
                   filename = "DEG_heatmap_wilcox.pdf")
  
  
  
  # Make heatmap for MAST analysis
  data_frame_list_MAST_DEG = list(L2_3_IT_MAST_DEG,
                                  L6_MAST_DEG,
                                  Sst_MAST_DEG,
                                  L5_MAST_DEG,
                                  L4_MAST_DEG,
                                  Pvalb_MAST_DEG,
                                  Sncg_MAST_DEG,
                                  Non_neuronal_MAST_DEG,
                                  Oligo_MAST_DEG,
                                  Vip_MAST_DEG,
                                  Lamp5_MAST_DEG,
                                  Astro_MAST_DEG,
                                  Peri_MAST_DEG,
                                  Endo_MAST_DEG)






################################################################################
# Limma Analysis
# By Osman Sharifi

clusterL2_3_IT <- subset(experiment.aggregate, idents = "L2_3_IT")
expr_L2_3_IT <- as.matrix(GetAssayData(clusterL2_3_IT))
# Filter out genes that are 0 for every cell in this cluster
bad_L2_3_IT <- which(rowSums(expr_L2_3_IT) == 0)
expr_L2_3_IT <- expr_L2_3_IT[-bad_L2_3_IT,]
mm_L2_3_IT <- model.matrix(~0 + orig.ident, data = clusterL2_3_IT@meta.data)
fitL2_3_IT <- lmFit(expr_L2_3_IT, mm_L2_3_IT)
head(coef(fitL2_3_IT)) # Means in each sample for each gene
contr_L2_3_IT<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitL2_3_IT)))
tmp_L2_3_IT <- contrasts.fit(fitL2_3_IT, contrasts = contr_L2_3_IT)
tmp_L2_3_IT <- eBayes(tmp_L2_3_IT)
L2_3_IT_toptable <- topTable(tmp_L2_3_IT, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(L2_3_IT_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_limma_top_20.csv")
L2_3_IT_Limma_stat_sig <- subset(x = L2_3_IT_toptable, subset = adj.P.Val < 0.05)
write.csv(L2_3_IT_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_Limma_DEG_only_stat_sig.csv")

clusterL6 <- subset(experiment.aggregate, idents = "L6")
expr_L6 <- as.matrix(GetAssayData(clusterL6))
# Filter out genes that are 0 for every cell in this cluster
bad_L6 <- which(rowSums(expr_L6) == 0)
expr_L6 <- expr_L6[-bad_L6,]
mm_L6 <- model.matrix(~0 + orig.ident, data = clusterL6@meta.data)
fitL6 <- lmFit(expr_L6, mm_L6)
head(coef(fitL6)) # Means in each sample for each gene
contr_L6<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitL6)))
tmp_L6 <- contrasts.fit(fitL6, contrasts = contr_L6)
tmp_L6 <- eBayes(tmp_L6)
L6_toptable <- topTable(tmp_L6, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(L6_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_limma_top_20.csv")
L6_Limma_stat_sig <- subset(x = L6_toptable, subset = adj.P.Val < 0.05)
write.csv(L6_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_Limma_DEG_only_stat_sig.csv")

clusterSst <- subset(experiment.aggregate, idents = "Sst")
expr_Sst <- as.matrix(GetAssayData(clusterSst))
# Filter out genes that are 0 for every cell in this cluster
bad_Sst <- which(rowSums(expr_Sst) == 0)
expr_Sst <- expr_Sst[-bad_Sst,]
mm_Sst <- model.matrix(~0 + orig.ident, data = clusterSst@meta.data)
fitSst <- lmFit(expr_Sst, mm_Sst)
head(coef(fitSst)) # Means in each sample for each gene
contr_Sst<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitSst)))
tmp_Sst <- contrasts.fit(fitSst, contrasts = contr_Sst)
tmp_Sst <- eBayes(tmp_Sst)
Sst_toptable <- topTable(tmp_Sst, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Sst_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_limma_top_20.csv")
Sst_Limma_stat_sig <- subset(x = Sst_toptable, subset = adj.P.Val < 0.05)
write.csv(Sst_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_Limma_DEG_only_stat_sig.csv")

clusterL5 <- subset(experiment.aggregate, idents = "L5")
expr_L5 <- as.matrix(GetAssayData(clusterL5))
# Filter out genes that are 0 for every cell in this cluster
bad_L5 <- which(rowSums(expr_L5) == 0)
expr_L5 <- expr_L5[-bad_L5,]
mm_L5 <- model.matrix(~0 + orig.ident, data = clusterL5@meta.data)
fitL5 <- lmFit(expr_L5, mm_L5)
head(coef(fitL5)) # Means in each sample for each gene
contr_L5<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitL5)))
tmp_L5 <- contrasts.fit(fitL5, contrasts = contr_L5)
tmp_L5 <- eBayes(tmp_L5)
L5_toptable <- topTable(tmp_L5, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(L5_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_limma_top_20.csv")
L5_Limma_stat_sig <- subset(x = L5_toptable, subset = adj.P.Val < 0.05)
write.csv(L5_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_Limma_DEG_only_stat_sig.csv")

clusterL4 <- subset(experiment.aggregate, idents = "L4")
expr_L4 <- as.matrix(GetAssayData(clusterL4))
# Filter out genes that are 0 for every cell in this cluster
bad_L4 <- which(rowSums(expr_L4) == 0)
expr_L4 <- expr_L4[-bad_L4,]
mm_L4 <- model.matrix(~0 + orig.ident, data = clusterL4@meta.data)
fitL4 <- lmFit(expr_L4, mm_L4)
head(coef(fitL4)) # Means in each sample for each gene
contr_L4<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitL4)))
tmp_L4 <- contrasts.fit(fitL4, contrasts = contr_L4)
tmp_L4 <- eBayes(tmp_L4)
L4_toptable <- topTable(tmp_L4, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(L4_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_limma_top_20.csv")
L4_Limma_stat_sig <- subset(x = L4_toptable, subset = adj.P.Val < 0.05)
write.csv(L4_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_Limma_DEG_only_stat_sig.csv")

clusterPvalb <- subset(experiment.aggregate, idents = "Pvalb")
expr_Pvalb <- as.matrix(GetAssayData(clusterPvalb))
# Filter out genes that are 0 for every cell in this cluster
bad_Pvalb <- which(rowSums(expr_Pvalb) == 0)
expr_Pvalb <- expr_Pvalb[-bad_Pvalb,]
mm_Pvalb <- model.matrix(~0 + orig.ident, data = clusterPvalb@meta.data)
fitPvalb <- lmFit(expr_Pvalb, mm_Pvalb)
head(coef(fitPvalb)) # Means in each sample for each gene
contr_Pvalb<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitPvalb)))
tmp_Pvalb <- contrasts.fit(fitPvalb, contrasts = contr_Pvalb)
tmp_Pvalb <- eBayes(tmp_Pvalb)
Pvalb_toptable <- topTable(tmp_Pvalb, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Pvalb_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_limma_top_20.csv")
Pvalb_Limma_stat_sig <- subset(x = Pvalb_toptable, subset = adj.P.Val < 0.05)
write.csv(Pvalb_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_Limma_DEG_only_stat_sig.csv")

clusterSncg <- subset(experiment.aggregate, idents = "Sncg")
expr_Sncg <- as.matrix(GetAssayData(clusterSncg))
# Filter out genes that are 0 for every cell in this cluster
bad_Sncg <- which(rowSums(expr_Sncg) == 0)
expr_Sncg <- expr_Sncg[-bad_Sncg,]
mm_Sncg <- model.matrix(~0 + orig.ident, data = clusterSncg@meta.data)
fitSncg <- lmFit(expr_Sncg, mm_Sncg)
head(coef(fitSncg)) # Means in each sample for each gene
contr_Sncg<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitSncg)))
tmp_Sncg <- contrasts.fit(fitSncg, contrasts = contr_Sncg)
tmp_Sncg <- eBayes(tmp_Sncg)
Sncg_toptable <- topTable(tmp_Sncg, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Sncg_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_limma_top_20.csv")
Sncg_Limma_stat_sig <- subset(x = Sncg_toptable, subset = adj.P.Val < 0.05)
write.csv(Sncg_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_Limma_DEG_only_stat_sig.csv")

clusterNon_neuronal <- subset(experiment.aggregate, idents = "Non-Neuronal")
expr_Non_neuronal <- as.matrix(GetAssayData(clusterNon_neuronal))
# Filter out genes that are 0 for every cell in this cluster
bad_Non_neuronal <- which(rowSums(expr_Non_neuronal) == 0)
expr_Non_neuronal <- expr_Non_neuronal[-bad_Non_neuronal,]
mm_Non_neuronal <- model.matrix(~0 + orig.ident, data = clusterNon_neuronal@meta.data)
fitNon_neuronal <- lmFit(expr_Non_neuronal, mm_Non_neuronal)
head(coef(fitNon_neuronal)) # Means in each sample for each gene
contr_Non_neuronal<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitNon_neuronal)))
tmp_Non_neuronal <- contrasts.fit(fitNon_neuronal, contrasts = contr_Non_neuronal)
tmp_Non_neuronal <- eBayes(tmp_Non_neuronal)
Non_neuronal_toptable <- topTable(tmp_Non_neuronal, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Non_neuronal_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_limma_top_20.csv")
Non_neuronal_Limma_stat_sig <- subset(x = Non_neuronal_toptable, subset = adj.P.Val < 0.05)
write.csv(Non_neuronal_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_Limma_DEG_only_stat_sig.csv")

clusterOligo <- subset(experiment.aggregate, idents = "Oligo")
expr_Oligo <- as.matrix(GetAssayData(clusterOligo))
# Filter out genes that are 0 for every cell in this cluster
bad_Oligo <- which(rowSums(expr_Oligo) == 0)
expr_Oligo <- expr_Oligo[-bad_Oligo,]
mm_Oligo <- model.matrix(~0 + orig.ident, data = clusterOligo@meta.data)
fitOligo <- lmFit(expr_Oligo, mm_Oligo)
head(coef(fitOligo)) # Means in each sample for each gene
contr_Oligo<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitOligo)))
tmp_Oligo <- contrasts.fit(fitOligo, contrasts = contr_Oligo)
tmp_Oligo <- eBayes(tmp_Oligo)
Oligo_toptable <- topTable(tmp_Oligo, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Oligo_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_limma_top_20.csv")
Oligo_Limma_stat_sig <- subset(x = Oligo_toptable, subset = adj.P.Val < 0.05)
write.csv(Oligo_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_Limma_DEG_only_stat_sig.csv")

clusterVip <- subset(experiment.aggregate, idents = "Vip")
expr_Vip <- as.matrix(GetAssayData(clusterVip))
# Filter out genes that are 0 for every cell in this cluster
bad_Vip <- which(rowSums(expr_Vip) == 0)
expr_Vip <- expr_Vip[-bad_Vip,]
mm_Vip <- model.matrix(~0 + orig.ident, data = clusterVip@meta.data)
fitVip <- lmFit(expr_Vip, mm_Vip)
head(coef(fitVip)) # Means in each sample for each gene
contr_Vip<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitVip)))
tmp_Vip <- contrasts.fit(fitVip, contrasts = contr_Vip)
tmp_Vip <- eBayes(tmp_Vip)
Vip_toptable <- topTable(tmp_Vip, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Vip_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_limma_top_20.csv")
Vip_Limma_stat_sig <- subset(x = Vip_toptable, subset = adj.P.Val < 0.05)
write.csv(Vip_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_Limma_DEG_only_stat_sig.csv")

clusterLamp5 <- subset(experiment.aggregate, idents = "Lamp5")
expr_Lamp5 <- as.matrix(GetAssayData(clusterLamp5))
# Filter out genes that are 0 for every cell in this cluster
bad_Lamp5 <- which(rowSums(expr_Lamp5) == 0)
expr_Lamp5 <- expr_Lamp5[-bad_Lamp5,]
mm_Lamp5 <- model.matrix(~0 + orig.ident, data = clusterLamp5@meta.data)
fitLamp5 <- lmFit(expr_Lamp5, mm_Lamp5)
head(coef(fitLamp5)) # Means in each sample for each gene
contr_Lamp5<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitLamp5)))
tmp_Lamp5 <- contrasts.fit(fitLamp5, contrasts = contr_Lamp5)
tmp_Lamp5 <- eBayes(tmp_Lamp5)
Lamp5_toptable <- topTable(tmp_Lamp5, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Lamp5_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_limma_top_20.csv")
Lamp5_Limma_stat_sig <- subset(x = Lamp5_toptable, subset = adj.P.Val < 0.05)
write.csv(Lamp5_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_Limma_DEG_only_stat_sig.csv")

clusterAstro <- subset(experiment.aggregate, idents = "Astro")
expr_Astro <- as.matrix(GetAssayData(clusterAstro))
# Filter out genes that are 0 for every cell in this cluster
bad_Astro <- which(rowSums(expr_Astro) == 0)
expr_Astro <- expr_Astro[-bad_Astro,]
mm_Astro <- model.matrix(~0 + orig.ident, data = clusterAstro@meta.data)
fitAstro <- lmFit(expr_Astro, mm_Astro)
head(coef(fitAstro)) # Means in each sample for each gene
contr_Astro<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitAstro)))
tmp_Astro <- contrasts.fit(fitAstro, contrasts = contr_Astro)
tmp_Astro <- eBayes(tmp_Astro)
Astro_toptable <- topTable(tmp_Astro, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Astro_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_limma_top_20.csv")
Astro_Limma_stat_sig <- subset(x = Astro_toptable, subset = adj.P.Val < 0.05)
write.csv(Astro_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_Limma_DEG_only_stat_sig.csv")

clusterPeri <- subset(experiment.aggregate, idents = "Peri")
expr_Peri <- as.matrix(GetAssayData(clusterPeri))
# Filter out genes that are 0 for every cell in this cluster
bad_Peri <- which(rowSums(expr_Peri) == 0)
expr_Peri <- expr_Peri[-bad_Peri,]
mm_Peri <- model.matrix(~0 + orig.ident, data = clusterPeri@meta.data)
fitPeri <- lmFit(expr_Peri, mm_Peri)
head(coef(fitPeri)) # Means in each sample for each gene
contr_Peri<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitPeri)))
tmp_Peri <- contrasts.fit(fitPeri, contrasts = contr_Peri)
tmp_Peri <- eBayes(tmp_Peri)
Peri_toptable <- topTable(tmp_Peri, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Peri_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_limma_top_20.csv")
Peri_Limma_stat_sig <- subset(x = Peri_toptable, subset = adj.P.Val < 0.05)
write.csv(Peri_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_Limma_DEG_only_stat_sig.csv")

clusterEndo <- subset(experiment.aggregate, idents = "Endo")
expr_Endo <- as.matrix(GetAssayData(clusterEndo))
# Filter out genes that are 0 for every cell in this cluster
bad_Endo <- which(rowSums(expr_Endo) == 0)
expr_Endo <- expr_Endo[-bad_Endo,]
mm_Endo <- model.matrix(~0 + orig.ident, data = clusterEndo@meta.data)
fitEndo <- lmFit(expr_Endo, mm_Endo)
head(coef(fitEndo)) # Means in each sample for each gene
contr_Endo<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitEndo)))
tmp_Endo <- contrasts.fit(fitEndo, contrasts = contr_Endo)
tmp_Endo <- eBayes(tmp_Endo)
Endo_toptable <- topTable(tmp_Endo, sort.by = "P", n = 20) # Top 20 DE genes
write.csv(Endo_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_limma_top_20.csv")
Endo_Limma_stat_sig <- subset(x = Endo_toptable, subset = adj.P.Val < 0.05)
write.csv(Endo_Limma_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_Limma_DEG_only_stat_sig.csv")

################################################################################
# DESeq2 Analysis
# By Viktoria Haghani

# Add one count to every RNA count so there are no zeroes in data set for DESeq2 log function (pseudocount)
# This is necessary because without pseudocounting, DESeq2 will have an error
experiment.aggregate[["RNA"]]@counts<-as.matrix(experiment.aggregate[["RNA"]]@counts)+1

L2_3_IT_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L2_3_IT", test.use = "DESeq2", slot = "counts")
write.csv("L2_3_IT_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_DESeq2_DEG_all_genes.csv")
L2_3_IT_DESeq2_DEG_stat_sig <- subset(x = L2_3_IT_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("L2_3_IT_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_DESeq2_DEG_stat_sig.csv")

L6_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L6", test.use = "DESeq2", slot = "counts")
write.csv("L6_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_DESeq2_DEG_all_genes.csv")
L6_DESeq2_DEG_stat_sig <- subset(x = L6_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("L6_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_DESeq2_DEG_stat_sig.csv")

Sst_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sst", test.use = "DESeq2", slot = "counts")
write.csv("Sst_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_DESeq2_DEG_all_genes.csv")
Sst_DESeq2_DEG_stat_sig <- subset(x = Sst_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Sst_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_DESeq2_DEG_stat_sig.csv")

L5_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L5", test.use = "DESeq2", slot = "counts")
write.csv("L5_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_DESeq2_DEG_all_genes.csv")
L5_DESeq2_DEG_stat_sig <- subset(x = L5_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("L5_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_DESeq2_DEG_stat_sig.csv")

L4_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L4", test.use = "DESeq2", slot = "counts")
write.csv("L4_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_DESeq2_DEG_all_genes.csv")
L4_DESeq2_DEG_stat_sig <- subset(x = L4_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("L4_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_DESeq2_DEG_stat_sig.csv")

Pvalb_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Pvalb", test.use = "DESeq2", slot = "counts")
write.csv("Pvalb_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_DESeq2_DEG_all_genes.csv")
Pvalb_DESeq2_DEG_stat_sig <- subset(x = Pvalb_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Pvalb_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_DESeq2_DEG_stat_sig.csv")

Sncg_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sncg", test.use = "DESeq2", slot = "counts")
write.csv("Sncg_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_DESeq2_DEG_all_genes.csv")
Sncg_DESeq2_DEG_stat_sig <- subset(x = Sncg_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Sncg_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_DESeq2_DEG_stat_sig.csv")

Non_neuronal_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Non_neuronal", test.use = "DESeq2", slot = "counts")
write.csv("Non_neuronal_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_DESeq2_DEG_all_genes.csv")
Non_neuronal_DESeq2_DEG_stat_sig <- subset(x = Non_neuronal_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Non_neuronal_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_DESeq2_DEG_stat_sig.csv")

Oligo_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Oligo", test.use = "DESeq2", slot = "counts")
write.csv("Oligo_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_DESeq2_DEG_all_genes.csv")
Oligo_DESeq2_DEG_stat_sig <- subset(x = Oligo_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Oligo_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_DESeq2_DEG_stat_sig.csv")

Vip_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Vip", test.use = "DESeq2", slot = "counts")
write.csv("Vip_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_DESeq2_DEG_all_genes.csv")
Vip_DESeq2_DEG_stat_sig <- subset(x = Vip_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Vip_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_DESeq2_DEG_stat_sig.csv")

Lamp5_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Lamp5", test.use = "DESeq2", slot = "counts")
write.csv("Lamp5_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_DESeq2_DEG_all_genes.csv")
Lamp5_DESeq2_DEG_stat_sig <- subset(x = Lamp5_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Lamp5_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_DESeq2_DEG_stat_sig.csv")

Astro_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Astro", test.use = "DESeq2", slot = "counts")
write.csv("Astro_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_DESeq2_DEG_all_genes.csv")
Astro_DESeq2_DEG_stat_sig <- subset(x = Astro_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Astro_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_DESeq2_DEG_stat_sig.csv")

Peri_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Peri", test.use = "DESeq2", slot = "counts")
write.csv("Peri_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_DESeq2_DEG_all_genes.csv")
Peri_DESeq2_DEG_stat_sig <- subset(x = Peri_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Peri_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_DESeq2_DEG_stat_sig.csv")

Endo_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Endo", test.use = "DESeq2", slot = "counts")
write.csv("Endo_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_DESeq2_DEG_all_genes.csv")
Endo_DESeq2_DEG_stat_sig <- subset(x = Endo_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Endo_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_DESeq2_DEG_stat_sig.csv")

################################################################################
# Venn Diagram for Differentially Expressed Genes Per Analysis
# By Viktoria Haghani

# List of genes differentially expressed per cluster for Limma
L2_3_IT_Limma_gene_list <- list(rownames(L2_3_IT_Limma_stat_sig))
L6_Limma_gene_list <- list(rownames(L6_Limma_stat_sig))
Sst_Limma_gene_list <- list(rownames(Sst_Limma_stat_sig))
L5_Limma_gene_list <- list(rownames(L5_Limma_stat_sig))
L4_Limma_gene_list <- list(rownames(L4_Limma_stat_sig))
Pvalb_Limma_gene_list <- list(rownames(Pvalb_Limma_stat_sig))
Sncg_Limma_gene_list <- list(rownames(Sncg_Limma_stat_sig))
Non_neuronal_Limma_gene_list <- list(rownames(Non_neuronal_Limma_stat_sig))
Oligo_Limma_gene_list <- list(rownames(Oligo_Limma_stat_sig))
Vip_Limma_gene_list <- list(rownames(Vip_Limma_stat_sig))
Lamp5_Limma_gene_list <- list(rownames(Lamp5_Limma_stat_sig))
Astro_Limma_gene_list <- list(rownames(Astro_Limma_stat_sig))
Peri_Limma_gene_list <- list(rownames(Peri_Limma_stat_sig))
Endo_Limma_gene_list <- list(rownames(Endo_Limma_stat_sig))


# List of genes differentially expressed per cluster for DESeq2
L2_3_IT_DESeq2_gene_list <- list()
L6_DESeq2_gene_list <- list()
Sst_DESeq2_gene_list <- list()
L5_DESeq2_gene_list <- list()
L4_DESeq2_gene_list <- list()
Pvalb_DESeq2_gene_list <- list()
Sncg_DESeq2_gene_list <- list()
Non_neuronal_DESeq2_gene_list <- list()
Oligo_DESeq2_gene_list <- list()
Vip_DESeq2_gene_list <- list()
Lamp5_DESeq2_gene_list <- list()
Astro_DESeq2_gene_list <- list()
Peri_DESeq2_gene_list <- list()
Endo_DESeq2_gene_list <- list()

# List of genes differentially expressed per cluster for EdgeR
L2_3_IT_EdgeR_gene_list <- list()
L6_EdgeR_gene_list <- list()
Sst_EdgeR_gene_list <- list()
L5_EdgeR_gene_list <- list()
L4_EdgeR_gene_list <- list()
Pvalb_EdgeR_gene_list <- list()
Sncg_EdgeR_gene_list <- list()
Non_neuronal_EdgeR_gene_list <- list()
Oligo_EdgeR_gene_list <- list()
Vip_EdgeR_gene_list <- list()
Lamp5_EdgeR_gene_list <- list()
Astro_EdgeR_gene_list <- list()
Peri_EdgeR_gene_list <- list()
Endo_EdgeR_gene_list <- list()

# Venn Diagram for Limma vs. DESeq2 per cluster
L2_3_IT_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L6_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sst_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L5_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L4_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Pvalb_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sncg_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Non_neuronal_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Oligo_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Vip_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Lamp5_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Astro_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Peri_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Endo_Limma_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")

# Venn Diagram for Limma vs. EdgeR per cluster
L2_3_IT_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L6_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sst_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L5_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L4_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Pvalb_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sncg_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Non_neuronal_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Oligo_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Vip_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Lamp5_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Astro_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Peri_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Endo_Limma_vs_EdgeR_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")

# Venn Diagram for EdgeR vs. DESeq2 per cluster
L2_3_IT_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L6_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sst_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L5_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L4_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Pvalb_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sncg_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Non_neuronal_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Oligo_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Vip_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Lamp5_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Astro_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Peri_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Endo_EdgeR_vs_DESeq2_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")

# Venn Diagram for Limma vs. DESeq2 vs. EdgeR per cluster
L2_3_IT_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L6_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sst_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L5_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
L4_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Pvalb_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Sncg_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Non_neuronal_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Oligo_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Vip_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Lamp5_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Astro_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Peri_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")
Endo_all_test_venn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram 1", nrtype="abs")