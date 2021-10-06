library(Seurat)
library(ggplot2)
library(topGO)
library(org.Mm.eg.db)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

## Paths
#data_file <- "/Users/osman/Desktop/LaSalle_lab/Scripts/P30_script/P30_Male_Cortex/P30_M_Cort_Labeled.RData"
data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P30_with_labels_proportions.rda"

## Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

################################################################################

## Load the Seurat object
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- "celltype.call"
options(width = 450)
table(Idents(experiment.aggregate),experiment.aggregate$orig.ident)

clusterL2_3_IT <- subset(experiment.aggregate, idents = 'L2_3_IT')
expr <- as.matrix(GetAssayData(clusterL2_3_IT))