library(Seurat)
library(edgeR)
library(limma)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
################################################################################
# EdgeR Analysis
# By Viktoria Haghani and Osman Sharifi
load('/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all_rett_mouse_cortex.RData')

# Subset the Seurat object
rett_mouse_cortex <- subset(all_rett_mouse_cortex, subset = Age != "E18")
DefaultAssay(rett_mouse_cortex) <- 'RNA'
Idents(rett_mouse_cortex) <- 'celltype.call'
#separate males and females
maleP30 <- subset(rett_mouse_cortex, subset = Age == "P30" & Sex == "Male")
maleP60 <- subset(rett_mouse_cortex, subset = Age == "P60" & Sex == "Male")
maleP120 <- subset(rett_mouse_cortex, subset = Age == "P120" & Sex == "Male")

femaleP30 <- subset(rett_mouse_cortex, subset = Age == "P30" & Sex == "Female")
femaleP60 <- subset(rett_mouse_cortex, subset = Age == "P60" & Sex == "Female")
femaleP120 <- subset(rett_mouse_cortex, subset = Age == "P150" & Sex == "Female")


setwd("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR_temporal")
# Assuming 'Condition' column indicates MUT and WT samples
# Replace 'Condition' with your actual column name if different

# Replace 'Condition' with your actual column name if different
condition <- maleP120$Condition

# List of cell types in your Seurat object
cell_types <- unique(Idents(maleP120))

# Loop over each cell type
for (cell_type in cell_types) {
  cat("Performing DEG analysis for cell type:", cell_type, "\n")
  
  # Subset Seurat object for the current cell type
  subset_obj <- subset(maleP120, idents = cell_type)
  
  # Get expression matrix
  counts <- as.matrix(subset_obj@assays$RNA@counts)
  
  # Get the metadata associated with the subset object
  subset_meta <- subset_obj@meta.data
  
  # Subsetting condition vector to include only samples for current cell type
  subset_condition <- condition[rownames(subset_meta)]
  
  # Create DGEList object
  dge <- DGEList(counts = counts, group = subset_condition)
  
  # Perform filtering and normalization
  dge <- calcNormFactors(dge)
  dge <- estimateCommonDisp(dge)
  
  # Perform differential expression analysis using classic mode
  dge <- glmFit(dge, design.matrix = model.matrix(~ subset_condition))
  dge <- glmLRT(dge, coef = 2)  # Coefficient for comparing MUT vs WT
  
  # Extract results
  results <- topTags(dge, n = Inf)
  
  # Save DEGs for MUT vs WT comparison
  degs <- results$table
  # You can customize this filename as needed
  write.csv(degs, file = paste0("DEGs_", cell_type, "_MUT_vs_WT.csv"))
  
  cat("DEG analysis for", cell_type, "complete.\n\n")
}
