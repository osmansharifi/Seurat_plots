library(Seurat)
library(limma)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths
data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P30_with_labels_proportions.rda"
DEG_data_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "Mice, Male, P30"

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
# Set threshold to 0.5%
experiment.aggregate <- subset(x = experiment.aggregate, subset = percent.mito <= "0.5")

################################################################################
## Limma Analysis

for (cell_type in cell_types){
}


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
L2_3_IT_toptable <- topTable(tmp_L2_3_IT, sort.by = "P", n = 50000) # Top 50000 DE genes (should cover all genes)
write.csv(L2_3_IT_toptable, file = "~/GitHub/snRNA-seq-pipeline/DEG_data/all_genes/L2_3_IT_limma_top_1000.csv")
L2_3_IT_Limma_stat_sig <- subset(x = L2_3_IT_toptable, subset = adj.P.Val < 0.05)
write.csv(L2_3_IT_Limma_stat_sig, file = '~/GitHub/snRNA-seq-pipeline/L2_3_IT_Limma_test.csv')

