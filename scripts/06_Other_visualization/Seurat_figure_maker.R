#ploting longitudinal data

library(Seurat)
library(scCustomize)
library(patchwork)
library(dplyr)

load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all.rett.combined.RData")
levels(all.rett.combined) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

#UMAP clustering plots
all_samples <- DimPlot_scCustom(all.rett.combined, label = FALSE)
all_samples_sex <- DimPlot_scCustom(all.rett.combined, group.by = "Sex", label = FALSE) 
all_samples_condition <- DimPlot_scCustom(all.rett.combined, group.by = "Condition", label = FALSE)
all_samples_age <- DimPlot_scCustom(all.rett.combined, group.by = "Age", label = FALSE)
all_samples_age+all_samples_condition+all_samples_sex
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/age_condition_sex_all_umap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
all_samples
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/allsamples_celltype.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# Cell proportions
cluster.averages <- AverageExpression(all.rett.combined)
head(cluster.averages[["RNA"]][,1:14])
cluster.averages <- AverageExpression(all.rett.combined, return.seurat = TRUE, group.by = "celltype.call")
CellScatter(cluster.averages, cell1 = "L2_3_IT", cell2 = "L4")
orig.levels <- levels(all.rett.combined)
Idents(all.rett.combined) <- gsub(pattern = " ", replacement = "_", x = Idents(all.rett.combined))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(all.rett.combined) <- orig.levels
cluster.averages <- AverageExpression(all.rett.combined, return.seurat = TRUE)
cluster.averages
# You can also plot heatmaps of these 'in silico' bulk datasets to visualize agreement between
# replicates
Idents(object = all.rett.combined) <- "orig.ident"
levels(all.rett.combined) <- c("WT_M_E18_WB1", "WT_M_E18_WB2", "MUT_M_E18_WB1", "MUT_M_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2", "MUT_F_E18_WB1", "MUT_F_E18_WB2","WT_M_P30_CORT1", "WT_M_P30_CORT2", "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_F_P30_CORT1", "WT_F_P30_CORT2", "MUT_F_P30_CORT1", "MUT_F_P30_CORT2","WT_M_P60_CORT1", "WT_M_P60_CORT2", "MUT_M_P60_CORT1", "MUT_M_P60_CORT2", "WT_F_P60_CORT1", "WT_F_P60_CORT2","MUT_F_P60_CORT1", "MUT_F_P60_CORT2", "WT_M_P120_CORT1", "WT_M_P120_CORT2", "MUT_M_P120_CORT1", "MUT_M_P120_CORT2", "WT_F_P150_CORT1", "WT_F_P150_CORT2", "WT_F_P150_CORT3", "WT_F_P150_CORT4",  "MUT_F_P150_CORT1", "MUT_F_P150_CORT2", "MUT_F_P150_CORT3", "MUT_F_P150_CORT4")

DoHeatmap(cluster.averages, features = unlist(TopFeatures(all.rett.combined[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette) + NoLegend()

DoHeatmap(cluster.averages, features = unlist(TopFeatures(all.rett.combined[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette, raster = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) 
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/replicate_heatmap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

DoHeatmap(cluster.averages, features = unlist(TopFeatures(all.rett.combined[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette) + scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

###Mecp2 expression
p1 <- FeaturePlot(object = all.rett.combined, features = "Mecp2", order = T)
p2 <- FeaturePlot_scCustom(seurat_object = all.rett.combined, features = "Mecp2", order = T)
p3 <- FeaturePlot_scCustom(seurat_object = all.rett.combined, colors_use = viridis_magma_dark_high, features = "Mecp2", order = T)
p4 <- FeaturePlot_scCustom(seurat_object = all.rett.combined, colors_use = viridis_light_high, features = "Mecp2", order = T)
p7 <- FeaturePlot_scCustom(seurat_object = all.rett.combined, features = 'WT_Mecp2', split.by = "Sex")
p8 <- FeaturePlot_scCustom(seurat_object = all.rett.combined, features = 'MUT_Mecp2', split.by = "Sex")
p5 <- FeaturePlot_scCustom(seurat_object = all.rett.combined, features = 'Mecp2', split.by = "Sex")
p6 <- FeaturePlot_scCustom(seurat_object = all.rett.combined, features = 'MUT_Mecp2', split.by = "Condition")
wrap_plots(p1, p2, p3, p4, ncol = 4) + plot_annotation(tag_levels = 'A')
wrap_plots(p5, p6, p7, p8, ncol = 2, nrow = 2) + plot_annotation(tag_levels = 'A')
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/Mecp2_expression.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
Idents(object = all.rett.combined) <- "celltype.call"
levels(all.rett.combined) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
cell_markers_manual <- c("Plch2","Sst","Vip", "Pvalb", "Slc17a8", "Macc1", "Rorb", "Fezf2", "Rprm", "Aqp4", "Rassf10", "Kcnj8", "Slc17a7", "Gad2", "Aspa")
Stacked_VlnPlot(all.rett.combined, features = cell_markers_manual, raster = FALSE, x_lab_rotate = TRUE)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/vlnplot_celltype.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
DoHeatmap(all.rett.combined, features = cell_markers_manual, size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette, raster = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) 
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/replicate_heatmap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

E18_subset <- subset(all.rett.combined, subset = orig.ident %in% c("WT_M_E18_WB1", "WT_M_E18_WB2", "MUT_M_E18_WB1", "MUT_M_E18_WB2"))

table(use.in$celltype)
