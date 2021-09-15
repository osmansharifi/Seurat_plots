library(ComplexHeatmap)
## Draws a heatmap with rows and columns seperated by a predefined label using ComplexHeatmap
plotDEG = function(data,       ## Data matrix
                   row.labels, ## Row annotation for data, assumed to be a vector of characters
                   col.labels, ## Column annotation for data, assumed to be a vector of characters
                   filename="~/P30_Male_cortex.png"){
  column_anno = HeatmapAnnotation(celltype=as.character(col.labels))
  ## Plot proportion heatmap
  DEG.heatmap = Heatmap(data,
                    top_annotation=column_anno,
                    col = c('purple', 'yellow'),
                    cluster_rows = T,
                    cluster_columns = T,
                    show_row_names = T,
                    row_split=factor(row.labels),
                    column_split=factor(col.labels),
                    show_heatmap_legend=T,
                    border=T)
  ## Plot
  png(filename, width=12, height=12, units='in', res=300)
  draw(DEG.heatmap)
  dev.off()
}
