library(Seurat)
library(edgeR)
library(limma)
library(ggplot2)
################################################################################
# EdgeR Analysis
# By Viktoria Haghani
load('/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all.rett.combined.RData')

# Subset the Seurat object
all.rett.combined_test <- subset(all.rett.combined, subset = Age == "E18")

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
