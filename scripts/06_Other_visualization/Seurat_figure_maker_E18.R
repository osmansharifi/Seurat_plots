#ploting longitudinal data

library(Seurat)
library(scCustomize)
library(patchwork)
library(dplyr)
library(Cairo)
library(Azimuth)
load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/mouse_embryonic_cortex.RData")
Idents(mouse_embryonic_cortex) <- "celltype.call"
levels(mouse_embryonic_cortex) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

#UMAP clustering plots
all_samples <- DimPlot_scCustom(mouse_embryonic_cortex, label = FALSE)
all_samples_sex <- DimPlot_scCustom(mouse_embryonic_cortex, group.by = "Sex", label = FALSE) 
all_samples_condition <- DimPlot_scCustom(mouse_embryonic_cortex, group.by = "Condition", label = FALSE)
all_samples_age <- DimPlot_scCustom(mouse_embryonic_cortex, group.by = "Age", label = FALSE)
all_samples+all_samples_age+all_samples_condition+all_samples_sex
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/age_condition_sex_all_umap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
all_samples+all_samples_age
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/allsamples_celltype_age.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
all_samples
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/allsamples_celltype.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
# Cell proportions
cluster.averages <- AverageExpression(mouse_embryonic_cortex)
head(cluster.averages[["RNA"]][,1:14])
cluster.averages <- AverageExpression(mouse_embryonic_cortex, return.seurat = TRUE, group.by = "celltype.call")
CellScatter(cluster.averages, cell1 = "L2_3_IT", cell2 = "L4")
orig.levels <- levels(mouse_embryonic_cortex)
Idents(mouse_embryonic_cortex) <- gsub(pattern = " ", replacement = "_", x = Idents(mouse_embryonic_cortex))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(mouse_embryonic_cortex) <- orig.levels
cluster.averages <- AverageExpression(mouse_embryonic_cortex, return.seurat = TRUE)
cluster.averages
# You can also plot heatmaps of these 'in silico' bulk datasets to visualize agreement between
# replicates
Idents(object = mouse_embryonic_cortex) <- "orig.ident"
levels(mouse_embryonic_cortex) <- c("WT_M_E18_WB1", "WT_M_E18_WB2", "MUT_M_E18_WB1", "MUT_M_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2", "MUT_F_E18_WB1", "MUT_F_E18_WB2","WT_M_P30_CORT1", "WT_M_P30_CORT2", "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_F_P30_CORT1", "WT_F_P30_CORT2", "MUT_F_P30_CORT1", "MUT_F_P30_CORT2","WT_M_P60_CORT1", "WT_M_P60_CORT2", "MUT_M_P60_CORT1", "MUT_M_P60_CORT2", "WT_F_P60_CORT1", "WT_F_P60_CORT2","MUT_F_P60_CORT1", "MUT_F_P60_CORT2", "WT_M_P120_CORT1", "WT_M_P120_CORT2", "MUT_M_P120_CORT1", "MUT_M_P120_CORT2", "WT_F_P150_CORT1", "WT_F_P150_CORT2", "WT_F_P150_CORT3", "WT_F_P150_CORT4",  "MUT_F_P150_CORT1", "MUT_F_P150_CORT2", "MUT_F_P150_CORT3", "MUT_F_P150_CORT4")

DoHeatmap(cluster.averages, features = unlist(TopFeatures(mouse_embryonic_cortex[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette) + NoLegend()

DoHeatmap(cluster.averages, features = unlist(TopFeatures(mouse_embryonic_cortex[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette, raster = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) 
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/replicate_heatmap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

DoHeatmap(cluster.averages, features = unlist(TopFeatures(mouse_embryonic_cortex[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette) + scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

###Mecp2 expression
p1 <- FeaturePlot(object = mouse_embryonic_cortex, features = "Mecp2", order = T)
p2 <- FeaturePlot_scCustom(seurat_object = mouse_embryonic_cortex, features = "Mecp2", order = T)
p3 <- FeaturePlot_scCustom(seurat_object = mouse_embryonic_cortex, colors_use = viridis_magma_dark_high, features = "Mecp2", order = T)
p4 <- FeaturePlot_scCustom(seurat_object = mouse_embryonic_cortex, colors_use = viridis_light_high, features = "Mecp2", order = T)
p7 <- FeaturePlot_scCustom(seurat_object = mouse_embryonic_cortex, features = 'WT_Mecp2', split.by = "Sex")
p8 <- FeaturePlot_scCustom(seurat_object = mouse_embryonic_cortex, features = 'MUT_Mecp2', split.by = "Sex")
p5 <- FeaturePlot_scCustom(seurat_object = mouse_embryonic_cortex, features = 'Mecp2', split.by = "Sex")
p6 <- FeaturePlot_scCustom(seurat_object = mouse_embryonic_cortex, features = 'MUT_Mecp2', split.by = "Condition")
wrap_plots(p1, p2, p3, p4, ncol = 4) + plot_annotation(tag_levels = 'A')
wrap_plots(p5, p6, p7, p8, ncol = 2, nrow = 2) + plot_annotation(tag_levels = 'A')
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/Mecp2_expression.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
Idents(object = mouse_embryonic_cortex) <- "celltype.call"
levels(mouse_embryonic_cortex) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
cell_markers_manual <- c("Plch2","Sst","Vip", "Pvalb", "Slc17a8", "Macc1", "Rorb", "Fezf2", "Rprm", "Aqp4", "Rassf10", "Kcnj8", "Slc17a7", "Gad2", "Aspa")
Stacked_VlnPlot(mouse_embryonic_cortex, features = cell_markers_manual, raster = FALSE, x_lab_rotate = TRUE)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/vlnplot_celltype.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
DoHeatmap(mouse_embryonic_cortex, features = cell_markers_manual, size = 3, draw.lines = FALSE, label = TRUE, group.colors = polychrome_palette, raster = FALSE) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) 
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/replicate_heatmap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
all_markers <- FindAllMarkers(object = test)
top5_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 5, named_vector = FALSE,
                                    make_unique = TRUE)

Clustered_DotPlot(seurat_object = test, features = top5_markers)
DotPlot_scCustom(seurat_object = test, features = top5_markers, flip_axes = TRUE, x_lab_rotate = TRUE) +
theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'plain', vjust = 1, family="Times", colour = 'black'),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="Times", colour = 'black'),
    plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1, family="Times", colour = 'black'),
    
    axis.text.x = element_text(angle = 60, size = 14, face = 'plain', hjust = 1.0, vjust = 1, family="Times", colour = 'black'),
    axis.text.y = element_text(angle = 0, size = 12, face = 'plain', vjust = 0.5, family="Times", colour = 'black'),
    axis.title = element_text(size = 14, face = 'plain', family="Times", colour = 'black'),
    axis.title.x = element_text(size = 14, face = 'plain', family="Times", colour = 'black', vjust = 1),
    axis.title.y = element_text(size = 14, face = 'plain', family="Times", colour = 'black'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold", family="Times"), # Text size
    title = element_text(size = 14, face = "plain", family="Times"))
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/celltype_markers_postnatal.pdf",
                device = NULL,
                height = 10,
                width = 12)



E18_subset <- subset(mouse_embryonic_cortex, subset = orig.ident %in% c("WT_M_E18_WB1", "WT_M_E18_WB2", "MUT_M_E18_WB1", "MUT_M_E18_WB2"))
postnat_subset <- subset(mouse_embryonic_cortex, subset = orig.ident %in% c("WT_M_P30_CORT1", "WT_M_P30_CORT2", "MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_F_P30_CORT1", "WT_F_P30_CORT2", "MUT_F_P30_CORT1", "MUT_F_P30_CORT2","WT_M_P60_CORT1", "WT_M_P60_CORT2", "MUT_M_P60_CORT1", "MUT_M_P60_CORT2", "WT_F_P60_CORT1", "WT_F_P60_CORT2","MUT_F_P60_CORT1", "MUT_F_P60_CORT2", "WT_M_P120_CORT1", "WT_M_P120_CORT2", "MUT_M_P120_CORT1", "MUT_M_P120_CORT2", "WT_F_P150_CORT1", "WT_F_P150_CORT2", "WT_F_P150_CORT3", "WT_F_P150_CORT4",  "MUT_F_P150_CORT1", "MUT_F_P150_CORT2", "MUT_F_P150_CORT3", "MUT_F_P150_CORT4"))
table(postnat_subset$celltype.call)
test <- RunUMAP(postnat_subset, dims = 1:20)
DimPlot_scCustom(postnat_subset, label = FALSE)
DimPlot_scCustom(test, label = FALSE)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/celltype_umap.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
levels(test) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
cell_markers_manual <- c("Plch2","Sst","Vip", "Pvalb", "Slc17a8", "Macc1", "Rorb", "Fezf2", "Rprm", "Aqp4", "Rassf10", "Kcnj8", "Slc17a7", "Gad2", "Aspa")
Stacked_VlnPlot(test, features = cell_markers_manual, raster = FALSE, x_lab_rotate = TRUE)
DotPlot_scCustom(test, features = cell_markers_manual, x_lab_rotate = TRUE, flip_axes = TRUE)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/dotplot_manual_markers.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
###Mecp2 expression
p1 <- FeaturePlot(object = test, features = "Mecp2", order = T)
p2 <- FeaturePlot_scCustom(seurat_object = test, features = "Mecp2", order = T)
p3 <- FeaturePlot_scCustom(seurat_object = test, colors_use = viridis_magma_dark_high, features = "Mecp2", order = T)
p4 <- FeaturePlot_scCustom(seurat_object = test, colors_use = viridis_light_high, features = "Mecp2", order = T)
p7 <- FeaturePlot_scCustom(seurat_object = test, features = 'WT_Mecp2', pt.size = 1)
p8 <- FeaturePlot_scCustom(seurat_object = test, features = 'MUT_Mecp2', pt.size = 1)
p5 <- FeaturePlot_scCustom(seurat_object = test, features = 'Mecp2', split.by = "Sex")
p6 <- FeaturePlot_scCustom(seurat_object = test, features = 'MUT_Mecp2', split.by = "Condition")
wrap_plots(p1, p2, p3, p4, ncol = 4) + plot_annotation(tag_levels = 'A')
wrap_plots(p7, p8, ncol = 2, nrow = 1) + plot_annotation(tag_levels = 'A')
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/postnat_Mecp2.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
