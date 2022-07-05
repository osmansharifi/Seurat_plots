##### Summary of scRNA-seq DEG Analysis #####

library(openxlsx)
library(glue)
library(tidyr)
library(ggplot2)
library(viridis)

s.obj.name = "rett_P30_with_labels_proportions"
setwd(glue::glue("/Users/karineier/Documents/scRNA-seq-differential-expression/{s.obj.name}"))

cell_types = list.dirs(glue::glue("/Users/karineier/Documents/scRNA-seq-differential-expression/{s.obj.name}/limmaVoomPB"))
cell_types = cell_types[-c(1, which(grepl("QC", cell_types)=="TRUE"), which(grepl("plotData", cell_types)=="TRUE"), which(grepl("interactivePlots", cell_types)=="TRUE"))]
#cell_types = cell_types[-c(which(grepl("", cell_types)=="TRUE"), 23)]


#cell_types.activated = cell_types[which(grepl("-activated", cell_types) == "TRUE")]
#cell_types.not = cell_types[which(grepl("-not", cell_types)=="TRUE")]
#cell_types.all = cell_types[-c(which(grepl("-activated", cell_types) == "TRUE"), which(grepl("-not", cell_types)=="TRUE"))]

list <- lapply(cell_types, function(cellType) {
  openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
})

cell_types = gsub(glue::glue("/Users/karineier/Documents/scRNA-seq-differential-expression/{s.obj.name}/limmaVoomPB/"), "", cell_types)

names(list) = cell_types

numDEGs = sapply(cell_types, function(cellType) {
  length(which(list[[cellType]][,4] < 0.1))
})

plotData = data.frame(cell_type = cell_types, numDEGs = as.numeric(numDEGs))

#plotData.activated = plotData[which(grepl("-activated", plotData$cell_type) == "TRUE"),]
#plotData.not = plotData[which(grepl("-not", plotData$cell_type)=="TRUE"),]
#plotData.activated.not = rbind(plotData.activated, plotData.not)
#plotData.combined = plotData[-c(which(grepl("-activated", plotData$cell_type) == "TRUE"), which(grepl("-not", plotData$cell_type)=="TRUE")),]

#plot.activated.not = ggplot(plotData.activated.not, aes(x=cell_type, y=numDEGs, fill=DEG_type)) +
#  geom_bar(stat="identity") +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  #scale_fill_viridis(option="plasma") +
#  xlab("") +
#  ylab("") 

#ggsave("Summary_Genotype_DEGs_by_celltype_females_P150_activated_unactivated_limmaVoomPB.pdf", width=8.5, height=11)


plot = ggplot(plotData, aes(x=cell_type, y=numDEGs, fill=numDEGs)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(legend.position="NULL", axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_viridis(option="plasma") +
  xlab("") +
  ylab("") 

ggsave(glue::glue("Summary_Genotype_DEGs_by_celltype_{s.obj.name}_limmaVoomPB.pdf"), width=8.5, height=11)




