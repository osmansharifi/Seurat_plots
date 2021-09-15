#Integrate all data frames together to create a DEG heatmap
data.frame.list = list(L2_3_IT_DEG, L4_DEG, L5_DEG, L6_DEG, Sst_DEG, Pvalb_DEG, Vip_DEG, Oligo_DEG, Lamp5_DEG, Nonneuronal_DEG)
#gene.set = unique(rownames(Reduce(rbind, data.frame.list)))
gene.set = unique(unlist(lapply(data.frame.list, rownames)))
ct.set   = unique(colnames(Reduce(cbind, data.frame.list)))
DEG.matrix = matrix(0, nrow=length(gene.set), ncol=length(ct.set), dimnames=list(gene.set, ct.set))
for(df in data.frame.list){
  DEG.matrix[rownames(df), colnames(df)] = df[,1] ## Not sure what here, pvalue? foldchange? or binary existence
}
DEG.df = as.data.frame(DEG.matrix)
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
