library(Seurat)
library(MAST)
library(limma)
library(edgeR)
library(DESeq2)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(ggVennDiagram)
library(edgeR)
library(udunits2)
library(scran)
library(Rcpp)

################################################################################
# Visualization and data preparation
# By Osman Sharifi & Viktoria Haghani

## Load in data
load('rett_P30_with_labels_proportions.rda')
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'
# Values represent cell numbers for each cell type
before_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident)
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

################################################################################
# EdgeR Analysis
# By Viktoria Haghani

# Subset Seurat object for the cluster
L2_3_IT_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "L2_3_IT")
# Generate a count matrix 
counts <- as.matrix(L2_3_IT_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- L2_3_IT_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
L2_3_IT_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(L2_3_IT_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
L6_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "L6")
# Generate a count matrix 
counts <- as.matrix(L6_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- L6_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
L6_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(L6_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Sst_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Sst")
# Generate a count matrix 
counts <- as.matrix(Sst_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Sst_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Sst_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Sst_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
L5_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "L5")
# Generate a count matrix 
counts <- as.matrix(L5_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- L5_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
L5_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(L5_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
L4_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "L4")
# Generate a count matrix 
counts <- as.matrix(L4_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- L4_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
L4_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(L4_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Pvalb_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Pvalb")
# Generate a count matrix 
counts <- as.matrix(Pvalb_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Pvalb_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Pvalb_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Pvalb_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Sncg_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Sncg")
# Generate a count matrix 
counts <- as.matrix(Sncg_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Sncg_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Sncg_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Sncg_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Non_neuronal_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Non_neuronal")
# Generate a count matrix 
counts <- as.matrix(Non_neuronal_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Non_neuronal_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Non_neuronal_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Non_neuronal_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Oligo_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Oligo")
# Generate a count matrix 
counts <- as.matrix(Oligo_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Oligo_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Oligo_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Oligo_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Vip_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Vip")
# Generate a count matrix 
counts <- as.matrix(Vip_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Vip_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Vip_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Vip_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Lamp5_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Lamp5")
# Generate a count matrix 
counts <- as.matrix(Lamp5_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Lamp5_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Lamp5_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Lamp5_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Astro_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Astro")
# Generate a count matrix 
counts <- as.matrix(Astro_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Astro_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Astro_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Astro_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Peri_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Peri")
# Generate a count matrix 
counts <- as.matrix(Peri_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Peri_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Peri_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Peri_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_EdgeR_DEG_only_stat_sig.csv")

# Subset Seurat object for the cluster
Endo_edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == "Endo")
# Generate a count matrix 
counts <- as.matrix(Endo_edgeR_obj@assays$RNA@counts)
counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
# Subset the meta data for filtered gene/cells
metadata <- Endo_edgeR_obj@meta.data
metadata <- metadata[,c("nCount_RNA", "orig.ident")]
metadata <- metadata[colnames(counts),]
# Make single cell experiment
sce <- SingleCellExperiment(assays = counts, 
                            colData = metadata)
# Convert to edgeR object
dge <- convertTo(sce, type="edgeR", assay.type = 1)
meta_dge <- dge$samples
meta_dge <- meta_dge[,c("lib.size","norm.factors")]
meta_dge <- cbind(meta_dge, metadata)
meta_dge$group <- factor(meta_dge$orig.ident)
levels(meta_dge$group)
meta_dge$group <- relevel(meta_dge$group, "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P30_CORT1", "WT_M_P30_CORT2")
dge$samples <- meta_dge
# Model fit
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group, data=dge$samples)
head(design)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
# Differential expression testing
my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P30_CORT1+groupWT_M_P30_CORT2) - c(groupMUT_M_P30_CORT1+groupMUT_M_P30_CORT2), levels = design)
qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
## All genes
# Use Benjamini-Hochberg correction for p-values
qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
# Note that topTags() outputs adjusted p-values in the FDR column
head(qlf.contrast.all.genes$table)
#write.csv(qlf.contrast.all.genes$table, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_EdgeR_DEG_all_genes.csv")
## Only statistically significant genes
Endo_EdgeR_stat_sig <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
write.csv(Endo_EdgeR_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_EdgeR_DEG_only_stat_sig.csv")

# Read in data to work with for EdgeR analysis (so analysis does not need to be rerun)
L2_3_IT_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/L2_3_IT_EdgeR_DEG_only_stat_sig.csv")
L6_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/L6_EdgeR_DEG_only_stat_sig.csv")
Sst_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Sst_EdgeR_DEG_only_stat_sig.csv")
L5_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/L5_EdgeR_DEG_only_stat_sig.csv")
L4_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/L4_EdgeR_DEG_only_stat_sig.csv")
Pvalb_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Pvalb_EdgeR_DEG_only_stat_sig.csv")
Sncg_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Sncg_EdgeR_DEG_only_stat_sig.csv")
Non_neuronal_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Non_neuronal_EdgeR_DEG_only_stat_sig.csv")
Oligo_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Oligo_EdgeR_DEG_only_stat_sig.csv")
Vip_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Vip_EdgeR_DEG_only_stat_sig.csv")
Lamp5_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Lamp5_EdgeR_DEG_only_stat_sig.csv")
Astro_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Astro_EdgeR_DEG_only_stat_sig.csv")
Peri_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Peri_EdgeR_DEG_only_stat_sig.csv")
Endo_EdgeR_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/Endo_EdgeR_DEG_only_stat_sig.csv")

################################################################################
# DESeq2 Analysis
# By Viktoria Haghani

# Add one count to every RNA count so there are no zeroes in data set for DESeq2 log function (pseudocount)
# This is necessary because without pseudocounting, DESeq2 will have an error
experiment.aggregate[["RNA"]]@counts<-as.matrix(experiment.aggregate[["RNA"]]@counts)+1

L2_3_IT_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L2_3_IT", test.use = "DESeq2", slot = "counts")
write.csv(L2_3_IT_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_DESeq2_DEG_all_genes.csv")
L2_3_IT_DESeq2_DEG_stat_sig <- subset(x = L2_3_IT_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(L2_3_IT_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L2_3_IT_DESeq2_DEG_stat_sig.csv")

L6_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L6", test.use = "DESeq2", slot = "counts")
write.csv(L6_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L6_DESeq2_DEG_all_genes.csv")
L6_DESeq2_DEG_stat_sig <- subset(x = L6_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(L6_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L6_DESeq2_DEG_stat_sig.csv")

Sst_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sst", test.use = "DESeq2", slot = "counts")
write.csv(Sst_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sst_DESeq2_DEG_all_genes.csv")
Sst_DESeq2_DEG_stat_sig <- subset(x = Sst_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Sst_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sst_DESeq2_DEG_stat_sig.csv")

L5_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L5", test.use = "DESeq2", slot = "counts")
write.csv(L5_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L5_DESeq2_DEG_all_genes.csv")
L5_DESeq2_DEG_stat_sig <- subset(x = L5_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(L5_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L5_DESeq2_DEG_stat_sig.csv")

L4_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "L4", test.use = "DESeq2", slot = "counts")
write.csv(L4_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L4_DESeq2_DEG_all_genes.csv")
L4_DESeq2_DEG_stat_sig <- subset(x = L4_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(L4_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/L4_DESeq2_DEG_stat_sig.csv")

Pvalb_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Pvalb", test.use = "DESeq2", slot = "counts")
write.csv(Pvalb_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Pvalb_DESeq2_DEG_all_genes.csv")
Pvalb_DESeq2_DEG_stat_sig <- subset(x = Pvalb_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Pvalb_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Pvalb_DESeq2_DEG_stat_sig.csv")

Sncg_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Sncg", test.use = "DESeq2", slot = "counts")
write.csv(Sncg_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Sncg_DESeq2_DEG_all_genes.csv")
Sncg_DESeq2_DEG_stat_sig <- subset(x = Sncg_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Sncg_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Sncg_DESeq2_DEG_stat_sig.csv")

Non_neuronal_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Non_neuronal", test.use = "DESeq2", slot = "counts")
write.csv(Non_neuronal_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Non_neuronal_DESeq2_DEG_all_genes.csv")
Non_neuronal_DESeq2_DEG_stat_sig <- subset(x = Non_neuronal_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Non_neuronal_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Non_neuronal_DESeq2_DEG_stat_sig.csv")

Oligo_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Oligo", test.use = "DESeq2", slot = "counts")
write.csv(Oligo_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Oligo_DESeq2_DEG_all_genes.csv")
Oligo_DESeq2_DEG_stat_sig <- subset(x = Oligo_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Oligo_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Oligo_DESeq2_DEG_stat_sig.csv")

Vip_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Vip", test.use = "DESeq2", slot = "counts")
write.csv(Vip_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Vip_DESeq2_DEG_all_genes.csv")
Vip_DESeq2_DEG_stat_sig <- subset(x = Vip_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Vip_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Vip_DESeq2_DEG_stat_sig.csv")

Lamp5_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Lamp5", test.use = "DESeq2", slot = "counts")
write.csv(Lamp5_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Lamp5_DESeq2_DEG_all_genes.csv")
Lamp5_DESeq2_DEG_stat_sig <- subset(x = Lamp5_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Lamp5_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Lamp5_DESeq2_DEG_stat_sig.csv")

Astro_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Astro", test.use = "DESeq2", slot = "counts")
write.csv(Astro_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Astro_DESeq2_DEG_all_genes.csv")
Astro_DESeq2_DEG_stat_sig <- subset(x = Astro_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Astro_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Astro_DESeq2_DEG_stat_sig.csv")

Peri_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Peri", test.use = "DESeq2", slot = "counts")
write.csv("Peri_DESeq2_DEG", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Peri_DESeq2_DEG_all_genes.csv")
Peri_DESeq2_DEG_stat_sig <- subset(x = Peri_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv("Peri_DESeq2_DEG_stat_sig", file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Peri_DESeq2_DEG_stat_sig.csv")

Endo_DESeq2_DEG <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30_CORT", group.by = "new.ident", subset.ident = "Endo", test.use = "DESeq2", slot = "counts")
write.csv(Endo_DESeq2_DEG, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/Endo_DESeq2_DEG_all_genes.csv")
Endo_DESeq2_DEG_stat_sig <- subset(x = Endo_DESeq2_DEG, subset = p_val_adj < 0.05)
write.csv(Endo_DESeq2_DEG_stat_sig, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/stat_sig/Endo_DESeq2_DEG_stat_sig.csv")

# Read in data to work with for DESeq2 analysis (so analysis does not need to be rerun)
L2_3_IT_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/L2_3_IT_DESeq2_DEG_stat_sig.csv")
L6_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/L6_DESeq2_DEG_stat_sig.csv")
Sst_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Sst_DESeq2_DEG_stat_sig.csv")
L5_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/L5_DESeq2_DEG_stat_sig.csv")
L4_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/L4_DESeq2_DEG_stat_sig.csv")
Pvalb_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Pvalb_DESeq2_DEG_stat_sig.csv")
Sncg_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Sncg_DESeq2_DEG_stat_sig.csv")
Non_neuronal_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Non_neuronal_DESeq2_DEG_stat_sig.csv")
Oligo_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Oligo_DESeq2_DEG_stat_sig.csv")
Vip_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Vip_DESeq2_DEG_stat_sig.csv")
Lamp5_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Lamp5_DESeq2_DEG_stat_sig.csv")
Astro_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Astro_DESeq2_DEG_stat_sig.csv")
Peri_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Peri_DESeq2_DEG_stat_sig.csv")
Endo_DESeq2_DEG_stat_sig <- read.csv(file = "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/Endo_DESeq2_DEG_stat_sig.csv")

################################################################################
# Venn Diagram for Differentially Expressed Genes Per Analysis
# By Viktoria Haghani

# Notes
# Sncg is excluded because no DEG were identified by any of the programs
# Oligo only has Limma and EdgeR since DESeq2 didn't identify any DEG
# Astro, Peri, and Endo only had DEG identified by one program each, so they are also excluded

# List of genes differentially expressed per cluster for Limma
L2_3_IT_Limma_gene_list <- L2_3_IT_Limma_stat_sig$X
L6_Limma_gene_list <- L6_Limma_stat_sig$X
Sst_Limma_gene_list <- Sst_Limma_stat_sig$X
L5_Limma_gene_list <- L5_Limma_stat_sig$X
L4_Limma_gene_list <- L4_Limma_stat_sig$X
Pvalb_Limma_gene_list <- Pvalb_Limma_stat_sig$X
Non_neuronal_Limma_gene_list <- Non_neuronal_Limma_stat_sig$X
Oligo_Limma_gene_list <- Oligo_Limma_stat_sig$X
Vip_Limma_gene_list <- Vip_Limma_stat_sig$X
Lamp5_Limma_gene_list <- Lamp5_Limma_stat_sig$X
Astro_Limma_gene_list <- Astro_Limma_stat_sig$X
Peri_Limma_gene_list <- Peri_Limma_stat_sig$X
Endo_Limma_gene_list <- Endo_Limma_stat_sig$X
unique_Limma_genes <- unique(c(L2_3_IT_Limma_stat_sig$X,
                     L6_Limma_stat_sig$X,
                     Sst_Limma_stat_sig$X,
                     L5_Limma_stat_sig$X, 
                     L4_Limma_stat_sig$X,
                     Pvalb_Limma_stat_sig$X,
                     Non_neuronal_Limma_stat_sig$X,
                     Oligo_Limma_stat_sig$X,
                     Vip_Limma_stat_sig$X,
                     Lamp5_Limma_stat_sig$X,
                     Astro_Limma_stat_sig$X,
                     Peri_Limma_stat_sig$X,
                     Endo_Limma_stat_sig$X))
all_Limma_genes_not_unique <- c(L2_3_IT_Limma_stat_sig$X,
                            L6_Limma_stat_sig$X,
                            Sst_Limma_stat_sig$X,
                            L5_Limma_stat_sig$X, 
                            L4_Limma_stat_sig$X,
                            Pvalb_Limma_stat_sig$X,
                            Non_neuronal_Limma_stat_sig$X,
                            Oligo_Limma_stat_sig$X,
                            Vip_Limma_stat_sig$X,
                            Lamp5_Limma_stat_sig$X,
                            Astro_Limma_stat_sig$X,
                            Peri_Limma_stat_sig$X,
                            Endo_Limma_stat_sig$X)

# List of genes differentially expressed per cluster for DESeq2
L2_3_IT_DESeq2_gene_list <- L2_3_IT_DESeq2_DEG_stat_sig$X
L6_DESeq2_gene_list <- L6_DESeq2_DEG_stat_sig$X
Sst_DESeq2_gene_list <- Sst_DESeq2_DEG_stat_sig$X
L5_DESeq2_gene_list <- L5_DESeq2_DEG_stat_sig$X
L4_DESeq2_gene_list <- L4_DESeq2_DEG_stat_sig$X
Pvalb_DESeq2_gene_list <- Pvalb_DESeq2_DEG_stat_sig$X
Non_neuronal_DESeq2_gene_list <- Non_neuronal_DESeq2_DEG_stat_sig$X
Oligo_DESeq2_gene_list <- Oligo_DESeq2_DEG_stat_sig$X
Vip_DESeq2_gene_list <- Vip_DESeq2_DEG_stat_sig$X
Lamp5_DESeq2_gene_list <- Lamp5_DESeq2_DEG_stat_sig$X
Astro_DESeq2_gene_list <- Astro_DESeq2_DEG_stat_sig$X
Peri_DESeq2_gene_list <- Peri_DESeq2_DEG_stat_sig$X
Endo_DESeq2_gene_list <- Endo_DESeq2_DEG_stat_sig$X
unique_DESeq2_genes <- unique(c(L2_3_IT_DESeq2_DEG_stat_sig$X,
                             L6_DESeq2_DEG_stat_sig$X,
                             Sst_DESeq2_DEG_stat_sig$X,
                             L5_DESeq2_DEG_stat_sig$X, 
                             L4_DESeq2_DEG_stat_sig$X,
                             Pvalb_DESeq2_DEG_stat_sig$X,
                             Non_neuronal_DESeq2_DEG_stat_sig$X,
                             Oligo_DESeq2_DEG_stat_sig$X,
                             Vip_DESeq2_DEG_stat_sig$X,
                             Lamp5_DESeq2_DEG_stat_sig$X,
                             Astro_DESeq2_DEG_stat_sig$X,
                             Peri_DESeq2_DEG_stat_sig$X,
                             Endo_DESeq2_DEG_stat_sig$X))

# List of genes differentially expressed per cluster for EdgeR
L2_3_IT_EdgeR_gene_list <- L2_3_IT_EdgeR_stat_sig$X
L6_EdgeR_gene_list <- L6_EdgeR_stat_sig$X
Sst_EdgeR_gene_list <- Sst_EdgeR_stat_sig$X
L5_EdgeR_gene_list <- L5_EdgeR_stat_sig$X
L4_EdgeR_gene_list <- L4_EdgeR_stat_sig$X
Pvalb_EdgeR_gene_list <- Pvalb_EdgeR_stat_sig$X
Non_neuronal_EdgeR_gene_list <- Non_neuronal_EdgeR_stat_sig$X
Oligo_EdgeR_gene_list <- Oligo_EdgeR_stat_sig$X
Vip_EdgeR_gene_list <- Vip_EdgeR_stat_sig$X
Lamp5_EdgeR_gene_list <- Lamp5_EdgeR_stat_sig$X
Astro_EdgeR_gene_list <- Astro_EdgeR_stat_sig$X
Peri_EdgeR_gene_list <- Peri_EdgeR_stat_sig$X
Endo_EdgeR_gene_list <- Endo_EdgeR_stat_sig$X
unique_EdgeR_genes <- unique(c(L2_3_IT_EdgeR_stat_sig$X,
                            L6_EdgeR_stat_sig$X,
                            Sst_EdgeR_stat_sig$X,
                            L5_EdgeR_stat_sig$X, 
                            L4_EdgeR_stat_sig$X,
                            Pvalb_EdgeR_stat_sig$X,
                            Non_neuronal_EdgeR_stat_sig$X,
                            Oligo_EdgeR_stat_sig$X,
                            Vip_EdgeR_stat_sig$X,
                            Lamp5_EdgeR_stat_sig$X,
                            Astro_EdgeR_stat_sig$X,
                            Peri_EdgeR_stat_sig$X,
                            Endo_EdgeR_stat_sig$X
                            ))

# Venn Diagram for Limma vs. DESeq2 vs. EdgeR per cluster
L2_3_IT_venn_list <- list(L2_3_IT_Limma_gene_list, L2_3_IT_EdgeR_gene_list, L2_3_IT_DESeq2_gene_list)
#L2_3_IT_venn <- 
ggVennDiagram(L2_3_IT_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L2_3_IT") +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("L2_3_IT_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

L6_venn_list <- list(L6_Limma_gene_list, L6_EdgeR_gene_list, L6_DESeq2_gene_list)
L6_venn <- ggVennDiagram(L6_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L6") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L6_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Sst_venn_list <- list(Sst_Limma_gene_list, Sst_EdgeR_gene_list, Sst_DESeq2_gene_list)
Sst_venn <- ggVennDiagram(Sst_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Sst") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Sst_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

L5_venn_list <- list(L5_Limma_gene_list, L5_EdgeR_gene_list, L5_DESeq2_gene_list)
L5_venn <- ggVennDiagram(L5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L5_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

L4_venn_list <- list(L4_Limma_gene_list, L4_EdgeR_gene_list, L4_DESeq2_gene_list)
L4_venn <- ggVennDiagram(L4_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L4") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L4_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Pvalb_venn_list <- list(Pvalb_Limma_gene_list, Pvalb_EdgeR_gene_list, Pvalb_DESeq2_gene_list)
Pvalb_venn <- ggVennDiagram(Pvalb_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Pvalb") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Pvalb_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Non_neuronal_venn_list <- list(Non_neuronal_Limma_gene_list, Non_neuronal_EdgeR_gene_list, Non_neuronal_DESeq2_gene_list)
Non_neuronal_venn <- ggVennDiagram(Non_neuronal_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Non_neuronal") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Non_neuronal_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Oligo_venn_list <- list(Oligo_Limma_gene_list, Oligo_EdgeR_gene_list)
Oligo_venn <- ggVennDiagram(Oligo_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Oligo") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Oligo_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Vip_venn_list <- list(Vip_Limma_gene_list, Vip_EdgeR_gene_list, Vip_DESeq2_gene_list)
Vip_venn <- ggVennDiagram(Vip_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Vip") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Vip_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Lamp5_venn_list <- list(Lamp5_Limma_gene_list, Lamp5_EdgeR_gene_list, Lamp5_DESeq2_gene_list)
Lamp5_venn <- ggVennDiagram(Lamp5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Lamp5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Lamp5_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

unique_venn_list <- list(unique_Limma_genes, unique_EdgeR_genes, unique_DESeq2_genes)
unique_genes_venn <- ggVennDiagram(unique_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Unique Differentially Expressed Genes Identified for All Cell Types") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("unique_genes_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

# Show genes identified by all methods for cell types
Reduce(intersect, list(L2_3_IT_Limma_gene_list, L2_3_IT_DESeq2_gene_list, L2_3_IT_EdgeR_gene_list))
Reduce(intersect, list(L6_Limma_gene_list, L6_DESeq2_gene_list, L6_EdgeR_gene_list))
Reduce(intersect, list(Sst_Limma_gene_list, Sst_DESeq2_gene_list, Sst_EdgeR_gene_list))
Reduce(intersect, list(L5_Limma_gene_list, L5_DESeq2_gene_list, L5_EdgeR_gene_list))
Reduce(intersect, list(L4_Limma_gene_list, L4_DESeq2_gene_list, L4_EdgeR_gene_list))
Reduce(intersect, list(Pvalb_Limma_gene_list, Pvalb_DESeq2_gene_list, Pvalb_EdgeR_gene_list))
Reduce(intersect, list(Vip_Limma_gene_list, Vip_DESeq2_gene_list, Vip_EdgeR_gene_list))
Reduce(intersect, list(Lamp5_Limma_gene_list, Lamp5_DESeq2_gene_list, Lamp5_EdgeR_gene_list))
