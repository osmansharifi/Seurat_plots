##### Summary of scRNA-seq DEG Analysis #####

library(openxlsx)
library(glue)
library(tidyr)
library(ggplot2)
library(viridis)

s.obj.name = "rett_P30_with_labels_proportions"
setwd(glue::glue("/Users/karineier/Documents/scRNA-seq-differential-expression/{s.obj.name}"))

cell_types = list.dirs(glue::glue("/Users/karineier/Documents/scRNA-seq-differential-expression/{s.obj.name}/limmaVoomCC"))
cell_types = cell_types[-c(1, which(grepl("QC", cell_types)=="TRUE"), which(grepl("plotData", cell_types)=="TRUE"), which(grepl("interactivePlots", cell_types)=="TRUE"))]
#cell_types = cell_types[-c(which(grepl("", cell_types)=="TRUE"), 23)]
cell_types

list <- lapply(cell_types, function(cellType) {
  openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
})

cell_types = gsub(glue::glue("/Users/karineier/Documents/scRNA-seq-differential-expression/{s.obj.name}/limmaVoomCC/"), "", cell_types)

names(list) = cell_types

numDEGs = sapply(cell_types, function(cellType) {
  length(which(list[[cellType]][,4] < 0.1))
})

plotData = data.frame(cell_type = cell_types, numDEGs = as.numeric(numDEGs))

plotData.activated = plotData[which(grepl("-activated", plotData$cell_type) == "TRUE"),]
plotData.not = plotData[which(grepl("-not", plotData$cell_type)=="TRUE"),]
plotData.activated.not = rbind(plotData.activated, plotData.not)
plotData.combined = plotData[-c(which(grepl("-activated", plotData$cell_type) == "TRUE"), which(grepl("-not", plotData$cell_type)=="TRUE")),]

plot.activated.not = ggplot(plotData.activated.not, aes(x=cell_type, y=numDEGs, fill=numDEGs)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_viridis(option="plasma") +
  xlab("") +
  ylab("") 

ggsave(glue::glue("Summary_Genotype_DEGs_by_celltype_{s.obj.name}_activated_unactivated_limmaVoomCC.pdf"), width=8.5, height=11)


plot.combined = ggplot(plotData.combined, aes(x=cell_type, y=numDEGs, fill=numDEGs)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(legend.position="NULL", axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_viridis(option="plasma") +
  xlab("") +
  ylab("") 

ggsave(glue::glue("Summary_Genotype_DEGs_by_celltype_{s.obj.name}_limmaVoomCC.pdf"), width=8.5, height=11)


## plot p-values by cell type

Pvalues = sapply(cell_types, function(x){
  list[[x]]$P.Value
})

Pvalues.vec = do.call(c, Pvalues)

name.lengths = sapply(cell_types, function(x){
  rep(glue::glue("{x}"), times = length(Pvalues[[x]]))
})
names.pvals = do.call(c, name.lengths)

data = data.frame(Pvalues = Pvalues.vec, cell_type=names.pvals)

data.activated = data[which(grepl("-activated", data$cell_type) == "TRUE"),]
data.not = data[which(grepl("-not", data$cell_type)=="TRUE"),]
data.activated.not = rbind(data.activated, data.not)
data.combined = data[-c(which(grepl("-activated", data$cell_type) == "TRUE"), which(grepl("-not", data$cell_type)=="TRUE")),]


ggplot(data=data.activated.not, aes(x=cell_type, y=-log2(Pvalues), color=cell_type)) +
  geom_point() +
  theme_classic() +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="None")

ggsave(glue::glue("Summary_Genotype_DEG_pvals_by_celltype_{s.obj.name}_activated_unactivated_limmaVoomCC.pdf"), width=8.5, height=11)

ggplot(data=data.combined, aes(x=cell_type, y=-log2(Pvalues), color=cell_type)) +
  geom_point() +
  theme_classic() +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="None")

ggsave(glue::glue("Summary_Genotype_DEG_pvals_by_celltype_{s.obj.name}_limmaVoomCC.pdf"), width=8.5, height=11)

