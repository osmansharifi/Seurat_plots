library(Seurat)
library(cowplot)
library(ggplot2)

# Visualization and data preparation
# By Osman Sharifi & Viktoria Haghani

################################################################################

data_file <- "~/GitHub/snRNA-seq-pipeline/data/rett_P30_with_labels_proportions.rda"
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

################################################################################

## Load in data
load('rett_P30_with_labels_proportions.rda')
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'
# Values represent cell numbers for each cell type
before_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident)
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')
# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')

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
# Call table to compare cells before G2M, S, and mt exclusion
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
# Call table to compare cells left after  G2M, S, and mt exclusion
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