library(scAlign)
library(SingleCellExperiment)

## Load in mouse
load("/share/lasallelab/Osman/scAlign_collab/processed/mouse.rda")

## Load in Osman data
load('/share/lasallelab/Osman/scAlign_collab/clusters_seurat_object.RData')
osman = experiment.aggregate.regress

## Common genes
mouse <- FindVariableFeatures(mouse, do.plot = F, nFeature=3000)
osman <- FindVariableFeatures(osman, do.plot = F, nFeature=3000)
genes.use = Reduce(intersect, list(VariableFeatures(mouse),
                                   VariableFeatures(osman),
                                   rownames(mouse),
                                   rownames(osman)))

## Create paired dataset SCE objects to pass into scAlignCreateObject
mouseSCE <- SingleCellExperiment(
    assays = list(counts = mouse@assays$RNA@counts[genes.use,],
                  logcounts  = mouse@assays$RNA@data[genes.use,],
                  scale.data = mouse@assays$RNA@scale.data[genes.use,])
)

osmanSCE <- SingleCellExperiment(
  assays = list(counts = osman@assays$RNA@counts[genes.use,],
                logcounts  = osman@assays$RNA@data[genes.use,],
                scale.data = osman@assays$RNA@scale.data[genes.use,])
)

scAlignRETT = scAlignCreateObject(sce.objects = list("mouse"=mouseSCE, "osman"=osmanSCE),
                                 #labels =
                                 data.use="scale.data",
                                 pca.reduce = TRUE,
                                 pcs.compute = 50,
                                 cca.reduce = TRUE,
                                 ccs.compute = 15,
                                 project.name = "scAlign_Osman_Rett")

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
