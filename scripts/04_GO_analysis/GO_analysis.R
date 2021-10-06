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
cell_types <- list("L2_3_IT")
#, "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 
topgo_ontologies <- list("BP", "CC", "MF")
################################################################################

## Load the Seurat object
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- "celltype.call"
options(width = 450)
# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')
table(Idents(experiment.aggregate),experiment.aggregate$orig.ident)

for (cell_type in cell_types){
  cell_cluster <- subset(experiment.aggregate, idents = cell_type)
  expr <- as.matrix(GetAssayData(cell_cluster))
  # Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
  n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
  expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.75)]
  all.genes <- rownames(expr)
  # Define geneList as 1 if gene is in expressed.genes, 0 otherwise
  geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
  names(geneList) <- all.genes
  for (ont in topgo_ontologies){
    # Create topGOdata object
    GOdata <- new("topGOdata",
                    ontology = ont,
                    allGenes = geneList,
                    geneSelectionFun = function(x)(x == 1),
                    annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
  }
}


