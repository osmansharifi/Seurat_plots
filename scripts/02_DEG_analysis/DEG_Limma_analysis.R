library(Seurat)
library(limma)
library(glue)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths

data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_E18_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P30_with_labels_proportions.rda"

DEG_data_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_M_E18_WB"

################################################################################
## Data preparation

# Load in data
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'

# Prepare data
# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')
# We want to get rid of the G2M and S phase cells, so subset to keep only G1 cells
experiment.aggregate <- subset(x = experiment.aggregate, subset = cell.cycle == "G1")
# Set mitochondrial threshold to 0.5%
experiment.aggregate <- subset(x = experiment.aggregate, subset = percent.mito <= "0.5")

################################################################################
## Limma Analysis

for (cell_type in cell_types){
  cluster_cell <- subset(experiment.aggregate, idents = cell_type)
  expr_cell <- as.matrix(GetAssayData(cluster_cell))
  # Filter out genes that are 0 for every cell in this cluster
  bad_cell <- which(rowSums(expr_cell) == 0)
  expr_cell <- expr_cell[-bad_cell,]
  mm_cell <- model.matrix(~0 + orig.ident, data = cluster_cell@meta.data)
  # Fit the model
  fit_cell <- lmFit(expr_cell, mm_cell)
  # Means in each sample for each gene
  head(coef(fit_cell)) 
  # Contrast WT-MUT accounting for repliicates
  contr_cell<- makeContrasts(c(orig.identWT_M_E18_WB1+orig.identWT_M_E18_WB2) - c(orig.identMUT_M_E18_WB1+orig.identMUT_M_E18_WB2), levels = colnames(coef(fit_cell)))
  tmp_cell <- contrasts.fit(fit_cell, contrasts = contr_cell)
  # Use empirical Bayes to calculate the t-statistics
  tmp_cell <- eBayes(tmp_cell)
  # Find top 50000 DE genes (should cover all genes)
  cell_toptable <- topTable(tmp_cell, sort.by = "P", n = 50000) 
  # Subset data to remove all non-significant genes
  cell_Limma_DEG <- subset(x = cell_toptable, subset = adj.P.Val < 0.05)
  # Write data to CSV so analysis does not need to be rerun when working with data
  write.csv(cell_Limma_DEG, file = glue(DEG_data_dir, cell_type, "_", metadata_info, "_Limma_DEG.csv"))
}
