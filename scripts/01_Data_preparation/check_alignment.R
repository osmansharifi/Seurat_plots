library(ggplot2)
library(cowplot)
library(Seurat)

## Load in mouse
load("/share/lasallelab/rett_female/scAlign_collab/processed/mouse.rda")
mouse <- obj
## Load in rett_female data
load('/share/lasallelab/rett_female/scAlign_collab/clusters_seurat_object.RData')

rett_female = experiment.aggregate

## Find var genes
mouse <- FindVariableFeatures(mouse, do.plot = F, nFeature=3000)
rett_female <- FindVariableFeatures(rett_female, do.plot = F, nFeature=3000)

## Combine our Seurat objects
combined.object <- merge(mouse, rett_female, add.cell.ids = c("ALLEN", "rett_female"), project = "RETT")
combined.object <- ScaleData(combined.object, do.scale=T, do.center=T, display.progress = T)

## Run PCA and UMAP
combined.object = RunPCA(combined.object)
combined.object = RunUMAP(combined.object, dims = 1:30)

## Plot tsne results
plot.me <- data.frame(x=combined.object@reductions$umap@cell.embeddings[,1],
                      y=combined.object@reductions$umap@cell.embeddings[,2],
                      labels=Idents(combined.object),
                      stringsAsFactors=FALSE)
unaligned.plot <- ggplot(plot.me, aes(x=x, y=y, colour = labels)) +
                  geom_point(size=2) +
                  scale_colour_manual(values=c("blue", "red")) +
                  xlab('UMAP_1') +
                  ylab('UMAP_2') +
                  theme_bw() +
                  theme(panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill = "transparent"), # bg of the panel
                        plot.background = element_rect(fill = "transparent", color = NA),
                        axis.line = element_line(colour = 'black',size=1))
plot(unaligned.plot)
ggplot2::ggsave("unaligned_allen_female.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
################################################################################
## Seurat alignment
################################################################################
rett.anchors <- FindIntegrationAnchors(object.list = list(mouse, rett_female), dims = 1:30)
rett.integrated <- IntegrateData(anchorset = rett.anchors, dims = 1:30)
DefaultAssay(rett.integrated) <- "integrated"

## Viz
rett.integrated <- ScaleData(rett.integrated, verbose = FALSE)
rett.integrated <- RunPCA(rett.integrated, npcs = 30, verbose = FALSE)
rett.integrated <- RunUMAP(rett.integrated, reduction = "pca", dims = 1:30)
## Plot
p1 <- DimPlot(rett.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(rett.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
## Save
png("/share/lasallelab/rett_female/scAlign_collab/seurat_alignment.png", width=12, height=12, res=150,  units="in")
plot_grid(p1, p2)
dev.off()

## Transfer labels
#https://satijalab.org/seurat/v3.1/integration.html
