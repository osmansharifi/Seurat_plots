library(Seurat)
library(patchwork)
library(ggplot2)

load("/Users/osman/Desktop/LaSalle_lab/Scripts/All_female_samples/all_female_cortex_labeled.RData")
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

features<- c("Gad1", "Slc17a7", "Vip", "Sst", "Pvalb", "Slc1a3", "Olig1", "Slc30a3", "Rorb")#, "Cux2")
StackedVlnPlot(obj = all_female.query, features = features, group.by = "predicted.id")

ggplot2::ggsave("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/seurat_figures/female/all_female_gene_marker_vln.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

mouse_colors_list <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid",
                       "orange", "gold", "gray", "red", "yellow", "black")
StackedVlnPlot(obj = all_female.query, features = features, colors_use = mouse_colors_list, flip = TRUE)
StackedVlnPlot(obj = all_female.query, features = features, group.by = "predicted.id", flip = TRUE)
