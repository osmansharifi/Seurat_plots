library(scAlign)
library(SingleCellExperiment)

## Load in mouse
load("/share/lasallelab/Osman/scAlign_collab/processed/mouse.rda")
#load("/Users/osman/Box Sync/single_nucleus_RNA-seq/rda_scAlign/mouse.rda")
mouse <- obj
## Load in rett_female data
load('/share/lasallelab/Osman/scAlign_collab/Cortex/all_female_preLabel.RData')
#load('/Users/osman/Desktop/LaSalle_lab/Scripts/All_female_samples/all_female_preLabel.RData')
rett_female = experiment.aggregate
rett_female.data = as.matrix(rett_female@assays$RNA@counts); rownames(rett_female.data) = toupper(rownames(rett_female.data))
rett_female = CreateSeuratObject(rett_female.data)
rett_female = NormalizeData(rett_female)
rett_female = ScaleData(rett_female)

## Common genes
mouse <- FindVariableFeatures(mouse, nFeature_RNA=3000)
rett_female <- FindVariableFeatures(rett_female, nFeature_RNA=3000)
genes.use = Reduce(intersect, list(VariableFeatures(mouse),
                                   VariableFeatures(rett_female),
                                   rownames(mouse),
                                   rownames(rett_female)))

## Create paired dataset SCE objects to pass into scAlignCreateObject
mouseSCE <- SingleCellExperiment(
    assays = list(counts = mouse@assays$RNA@counts[genes.use,],
                  logcounts  = mouse@assays$RNA@data[genes.use,],
                  scale.data = mouse@assays$RNA@scale.data[genes.use,])
)

rett_femaleSCE <- SingleCellExperiment(
  assays = list(counts = rett_female@assays$RNA@counts[genes.use,],
                logcounts  = rett_female@assays$RNA@data[genes.use,],
                scale.data = rett_female@assays$RNA@scale.data[genes.use,])
)

scAlignRETT = scAlignCreateObject(sce.objects = list("mouse"=mouseSCE, "rett_female"=rett_femaleSCE),
                                 #labels = list(labels[which("mouse" = mouseSCE)],
                                              # lables[which("rett_female" = rett_femaleSCE)]),
                                 data.use="scale.data",
                                 pca.reduce = TRUE,
                                 pcs.compute = 50,
                                 cca.reduce = TRUE,
                                 ccs.compute = 15,
                                 project.name = "scAlign_female_Rett")

## Run scAlign with high_var_genes as input to the encoder (alignment) and logcounts with the decoder (projections).
scAlignRETT = scAlign(scAlignRETT,
                      options=scAlignOptions(steps=5000, log.every=5000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture="small"),
                      encoder.data="scale.data",
                      #decoder.data="logcounts",
                      supervised='none',
                      run.encoder=TRUE,
                      #run.decoder=TRUE,
                      #log.dir=file.path(results.dir, 'models','gene_input'),
                      device="GPU")
