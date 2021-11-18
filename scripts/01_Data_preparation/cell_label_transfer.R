# load libraries
library(scAlign)

# run alignment function
if (FALSE) {
#load Seurat Objects
load("/share/lasallelab/Osman/test_alignment/all_male_samples.RData")
all_male.reference <- all_male
load("/share/lasallelab/Osman/2021_PEBBLES_Cortex/PCB_FEMALE_CLUSTERS.RData")

pbmc.query <- pbmc3k[, 1351:2700]

# perform standard preprocessing on each object
pbmc.reference <- NormalizeData(pbmc.reference)
pbmc.reference <- FindVariableFeatures(pbmc.reference)
pbmc.reference <- ScaleData(pbmc.reference)

pbmc.query <- NormalizeData(pbmc.query)
pbmc.query <- FindVariableFeatures(pbmc.query)
pbmc.query <- ScaleData(pbmc.query)

# find anchors
anchors <- FindTransferAnchors(reference = pbmc.reference, query = pbmc.query)

# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = pbmc.reference$seurat_annotations)
pbmc.query <- AddMetaData(object = pbmc.query, metadata = predictions)
}