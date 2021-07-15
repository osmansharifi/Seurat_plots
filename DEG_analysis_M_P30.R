library(Seurat)
library(DESeq2)
library(limma)
library(MAST)
library(glue)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)

################################################################################
# Visualization and data preparation
# By Osman Sharifi & Viktoria Haghani

## Load in data
load('/Users/osman/Desktop/LaSalle_lab/Scripts/P30_script/P30_Male_Cortex/rett_P30_with_labels_proportions.rda')
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'
# Values represent cell numbers for each cell type
before_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident) 

## Subset cells in G1 and visualize UMAP
# Visualize UMAP
pcs.use <- 10
experiment.aggregate<- RunUMAP(experiment.aggregate, dims = 1:pcs.use)
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
DimPlot(experiment.aggregate, reduction = "umap", group.by = "celltype.call", label = TRUE) +
  ggtitle("Cell Types After G1 & mt Subsetting")
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
    file_name <- glue('{cell_type}_{test}_DEG')
    DEG_data[[length(DEG_data) + 1]] <- file_name
    file <- glue('{file_name}_only_stat_sig.csv')
    # This writes directly to a value called "file_name" instead of the variable; need to fix it later
    file_name <- FindMarkers(experiment.aggregate, ident.1 = "MUT_M_P30", group.by = "new.ident", subset.ident = cell_type, test.use = test)
    #file_name <- subset(x = file_name, subset = p_val_adj < 0.05)
    write.csv(file_name, file = file)
  }}

# Read all data from csv files so analysis does not need to be rerun every time
# Remove genes that are not statistically significant
# Write a for-loop for this eventually
L2_3_IT_wilcox_DEG <- read.csv(file = 'L2_3_IT_wilcox_DEG.csv')
L2_3_IT_wilcox_DEG <- subset(x = L2_3_IT_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L2_3_IT_wilcox_DEG, file = "L2_3_IT_wilcox_DEG_only_stat_sig.csv")
L2_3_IT_wilcox_DEG_df <- data.frame(data_frame_list_wilcox_DEG[,c(0,1)])

L6_wilcox_DEG <- read.csv(file = 'L6_wilcox_DEG.csv')
L6_wilcox_DEG <- subset(x = L6_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L6_wilcox_DEG, file = "L6_wilcox_DEG_only_stat_sig.csv")

Sst_wilcox_DEG <- read.csv(file = 'Sst_wilcox_DEG.csv')
Sst_wilcox_DEG <- subset(x = Sst_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Sst_wilcox_DEG, file = "Sst_wilcox_DEG_only_stat_sig.csv")

L5_wilcox_DEG <- read.csv(file = 'L5_wilcox_DEG.csv')
L5_wilcox_DEG <- subset(x = L5_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L5_wilcox_DEG, file = "L5_wilcox_DEG_only_stat_sig.csv")

L4_wilcox_DEG <- read.csv(file = 'L4_wilcox_DEG.csv')
L4_wilcox_DEG <- subset(x = L4_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(L4_wilcox_DEG, file = "L4_wilcox_DEG_only_stat_sig.csv")

Pvalb_wilcox_DEG <- read.csv(file = 'Pvalb_wilcox_DEG.csv')
Pvalb_wilcox_DEG <- subset(x = Pvalb_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Pvalb_wilcox_DEG, file = "Pvalb_wilcox_DEG_only_stat_sig.csv")

Sncg_wilcox_DEG <- read.csv(file = 'Sncg_wilcox_DEG.csv')
Sncg_wilcox_DEG <- subset(x = Sncg_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Sncg_wilcox_DEG, file = "Sncg_wilcox_DEG_only_stat_sig.csv")

Non_neuronal_wilcox_DEG <- read.csv(file = 'Non_neuronal_wilcox_DEG.csv')
Non_neuronal_wilcox_DEG <- subset(x = Non_neuronal_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Non_neuronal_wilcox_DEG, file = "Non_neuronal_wilcox_DEG_only_stat_sig.csv")

Oligo_wilcox_DEG <- read.csv(file = 'Oligo_wilcox_DEG.csv')
Oligo_wilcox_DEG <- subset(x = Oligo_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Oligo_wilcox_DEG, file = "Oligo_wilcox_DEG_only_stat_sig.csv")

Vip_wilcox_DEG <- read.csv(file = 'Vip_wilcox_DEG.csv')
Vip_wilcox_DEG <- subset(x = Vip_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Vip_wilcox_DEG, file = "Vip_wilcox_DEG_only_stat_sig.csv")

Lamp5_wilcox_DEG <- read.csv(file = 'Lamp5_wilcox_DEG.csv')
Lamp5_wilcox_DEG <- subset(x = Lamp5_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Lamp5_wilcox_DEG, file = "Lamp5_wilcox_DEG_only_stat_sig.csv")

Astro_wilcox_DEG <- read.csv(file = 'Astro_wilcox_DEG.csv')
Astro_wilcox_DEG <- subset(x = Astro_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Astro_wilcox_DEG, file = "Astro_wilcox_DEG_only_stat_sig.csv")

Peri_wilcox_DEG <- read.csv(file = 'Peri_wilcox_DEG.csv')
Peri_wilcox_DEG <- subset(x = Peri_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Peri_wilcox_DEG, file = "Peri_wilcox_DEG_only_stat_sig.csv")

Endo_wilcox_DEG <- read.csv(file = 'Endo_wilcox_DEG.csv')
Endo_wilcox_DEG <- subset(x = Endo_wilcox_DEG, subset = p_val_adj < 0.05)
write.csv(Endo_wilcox_DEG, file = "Endo_wilcox_DEG_only_stat_sig.csv")

L2_3_IT_MAST_DEG <- read.csv(file = 'L2_3_IT_MAST_DEG.csv')
L2_3_IT_MAST_DEG <- subset(x = L2_3_IT_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L2_3_IT_MAST_DEG, file = "L2_3_IT_MAST_DEG_only_stat_sig.csv")

L6_MAST_DEG <- read.csv(file = 'L6_MAST_DEG.csv')
L6_MAST_DEG <- subset(x = L6_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L6_MAST_DEG, file = "L6_MAST_DEG_only_stat_sig.csv")

Sst_MAST_DEG <- read.csv(file = 'Sst_MAST_DEG.csv')
Sst_MAST_DEG <- subset(x = Sst_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Sst_MAST_DEG, file = "Sst_MAST_DEG_only_stat_sig.csv")

L5_MAST_DEG <- read.csv(file = 'L5_MAST_DEG.csv')
L5_MAST_DEG <- subset(x = L5_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L5_MAST_DEG, file = "L5_MAST_DEG_only_stat_sig.csv")

L4_MAST_DEG <- read.csv(file = 'L4_MAST_DEG.csv')
L4_MAST_DEG <- subset(x = L4_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(L4_MAST_DEG, file = "L4_MAST_DEG_only_stat_sig.csv")

Pvalb_MAST_DEG <- read.csv(file = 'Pvalb_MAST_DEG.csv')
Pvalb_MAST_DEG <- subset(x = Pvalb_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Pvalb_MAST_DEG, file = "Pvalb_MAST_DEG_only_stat_sig.csv")

Sncg_MAST_DEG <- read.csv(file = 'Sncg_MAST_DEG.csv')
Sncg_MAST_DEG <- subset(x = Sncg_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Sncg_MAST_DEG, file = "Sncg_MAST_DEG_only_stat_sig.csv")

Non_neuronal_MAST_DEG <- read.csv(file = 'Non_neuronal_MAST_DEG.csv')
Non_neuronal_MAST_DEG <- subset(x = Non_neuronal_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Non_neuronal_MAST_DEG, file = "Non_neuronal_MAST_DEG_only_stat_sig.csv")

Oligo_MAST_DEG <- read.csv(file = 'Oligo_MAST_DEG.csv')
Oligo_MAST_DEG <- subset(x = Oligo_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Oligo_MAST_DEG, file = "Oligo_MAST_DEG_only_stat_sig.csv")

Vip_MAST_DEG <- read.csv(file = 'Vip_MAST_DEG.csv')
Vip_MAST_DEG <- subset(x = Vip_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Vip_MAST_DEG, file = "Vip_MAST_DEG_only_stat_sig.csv")

Lamp5_MAST_DEG <- read.csv(file = 'Lamp5_MAST_DEG.csv')
Lamp5_MAST_DEG <- subset(x = Lamp5_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Lamp5_MAST_DEG, file = "Lamp5_MAST_DEG_only_stat_sig.csv")

Astro_MAST_DEG <- read.csv(file = 'Astro_MAST_DEG.csv')
Astro_MAST_DEG <- subset(x = Astro_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Astro_MAST_DEG, file = "Astro_MAST_DEG_only_stat_sig.csv")

Peri_MAST_DEG <- read.csv(file = 'Peri_MAST_DEG.csv')
Peri_MAST_DEG <- subset(x = Peri_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Peri_MAST_DEG, file = "Peri_MAST_DEG_only_stat_sig.csv")

Endo_MAST_DEG <- read.csv(file = 'Endo_MAST_DEG.csv')
Endo_MAST_DEG <- subset(x = Endo_MAST_DEG, subset = p_val_adj < 0.05)
write.csv(Endo_MAST_DEG, file = "Endo_MAST_DEG_only_stat_sig.csv")


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

clusterAstro <- subset(experiment.aggregate, idents = 'Astro')
expr_Astro <- as.matrix(GetAssayData(clusterAstro))

# Filter out genes that are 0 for every cell in this cluster
bad_Astro <- which(rowSums(expr_Astro) == 0)
expr_Astro <- expr_Astro[-bad_Astro,]

mm_Astro <- model.matrix(~0 + orig.ident, data = clusterAstro@meta.data)
fitAstro <- lmFit(expr_Astro, mm_Astro)  
head(coef(fitAstro)) # means in each sample for each gene
contr_Astro<- makeContrasts(c(orig.identWT_M_P30_CORT1+orig.identWT_M_P30_CORT2) - c(orig.identMUT_M_P30_CORT1+orig.identMUT_M_P30_CORT2), levels = colnames(coef(fitAstro)))
tmp_Astro <- contrasts.fit(fitAstro, contrasts = contr_Astro)
tmp_Astro <- eBayes(tmp_Astro)
topTable(tmp_Astro, sort.by = "P", n = 20) # top 20 DE genes
# HErere is a random change 

################################################################################
# DESeq2 Analysis
# By Viktoria Haghani

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