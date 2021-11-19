# load libraries
library(scAlign)

# run alignment function
if (FALSE) {
#load Seurat Objects
load("/share/lasallelab/Osman/test_alignment/all_male_samples.RData")
all_male.reference <- all_male
load("/share/lasallelab/Osman/2021_PEBBLES_Cortex/PCB_FEMALE_CLUSTERS.RData")
PCB_female.query <- experiment.aggregate

# perform standard preprocessing on each object
all_male.reference <- NormalizeData(all_male.reference)
all_male.reference <- FindVariableFeatures(all_male.reference)
all_male.reference <- ScaleData(all_male.reference)

PCB_female.query <- NormalizeData(PCB_female.query)
PCB_female.query <- FindVariableFeatures(PCB_female.query)
PCB_female.query <- ScaleData(PCB_female.query)

# find anchors
anchors <- FindTransferAnchors(reference = all_male.reference, query = PCB_female.query)

# transfer labels
celltype <- TransferData(anchorset = anchors, refdata = all_male.reference$celltype.call)
PCB_female.query <- AddMetaData(object = PCB_female.query, metadata = celltype)
}
all_female_pebbles <- PCB_female.query
save(all_female_pebbles, file="all_female_pebbles_labeled.RData")
