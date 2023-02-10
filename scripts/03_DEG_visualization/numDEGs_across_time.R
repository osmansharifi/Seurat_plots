### Summary Plots Across Time ###

library("ggplot2")
library("tidyr")
library("dplyr")

#Custom color palette
polychrome_palette  <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

# Order that celltypes should appear in 
x = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")

###############
## load male ##
###############
s.obj.names = c("M_MUT_and_WT_M_P30_CORT", "M_MUT_and_WT_M_P60_CORT", "M_MUT_and_WT_M_P120_CORT")
DEGlists = vector(mode="list", length=3)
plotDataList = vector(mode="list", length=3)

names(DEGlists) = s.obj.names
names(plotDataList) = s.obj.names

for(i in s.obj.names) {
  
  setwd(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/DEG_visualization"))
  
  cell_types = list.dirs(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/DEsingle"))
  cell_types = cell_types[-c(1, which(grepl("QC", cell_types)=="TRUE"), which(grepl("plotData", cell_types)=="TRUE"), which(grepl("interactivePlots", cell_types)=="TRUE"))]
  
  DEGlists[[i]] <- lapply(cell_types, function(cellType) {
    openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
  })
  
  cell_types = gsub(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/DEsingle/"), "", cell_types)
  
  names(DEGlists[[i]]) = cell_types
  
  numDEGs = sapply(cell_types, function(cellType) {
    length(which(DEGlists[[i]][[cellType]][,4] < 0.1))
  })
  
  plotDataList[[i]] = data.frame(cell_type = cell_types, numDEGs = as.numeric(numDEGs))
  
}

### Number of DEGs across time
#female
#plotDataList$M_MUT_and_WT_F_E18_CORT$Age = c(rep(-2, times=length(plotDataList$M_MUT_and_WT_F_E18_CORT$cell_type)))
#male
#plotDataList$M_MUT_and_WT_M_E18_CORT$Age = c(rep(-2, times=length(plotDataList$M_MUT_and_WT_M_E18_CORT$cell_type)))
plotDataList$M_MUT_and_WT_M_P30_CORT$Age = c(rep(30, times=length(plotDataList$M_MUT_and_WT_M_P30_CORT$cell_type)))
plotDataList$M_MUT_and_WT_M_P60_CORT$Age = c(rep(60, times=length(plotDataList$M_MUT_and_WT_M_P60_CORT$cell_type)))
plotDataList$M_MUT_and_WT_M_P120_CORT$Age = c(rep(120, times=length(plotDataList$M_MUT_and_WT_M_P120_CORT$cell_type)))
plotDataTogether = do.call(rbind, plotDataList)

plotData = plotDataTogether[-c(1, which(grepl("-activated", plotDataTogether$cell_type) == "TRUE"), which(grepl("-not", plotDataTogether$cell_type) == "TRUE")),]

plotData_male_desingle <- plotData %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type) 

male_desingle <- ggplot(data=plotData_male_desingle, aes(x=Age, y=numDEGs, color=cell_type)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = polychrome_palette)+
  geom_point(size = 3) +
  theme_classic() +
  xlab("Age (days)") +
  ggtitle("DEsingle Postnatal male numDEGs")+
  scale_x_continuous(breaks=c(30, 60, 120))+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'none',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 14, face = 'bold', hjust = 1.0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 14, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/numDEGs_male_postnatal_DEsingle.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

for(i in s.obj.names) {
  
  setwd(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/DEG_visualization"))
  
  cell_types = list.dirs(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/limmaVoomCC"))
  cell_types = cell_types[-c(1, which(grepl("QC", cell_types)=="TRUE"), which(grepl("plotData", cell_types)=="TRUE"), which(grepl("interactivePlots", cell_types)=="TRUE"))]
  
  DEGlists[[i]] <- lapply(cell_types, function(cellType) {
    openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
  })
  
  cell_types = gsub(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/limmaVoomCC/"), "", cell_types)
  
  names(DEGlists[[i]]) = cell_types
  
  numDEGs = sapply(cell_types, function(cellType) {
    length(which(DEGlists[[i]][[cellType]][,4] < 0.05))
  })
  
  plotDataList[[i]] = data.frame(cell_type = cell_types, numDEGs = as.numeric(numDEGs))
  
}

### Number of DEGs across time
#female
#plotDataList$M_MUT_and_WT_F_E18_CORT$Age = c(rep(-2, times=length(plotDataList$M_MUT_and_WT_F_E18_CORT$cell_type)))
#male
#plotDataList$M_MUT_and_WT_M_E18_CORT$Age = c(rep(-2, times=length(plotDataList$M_MUT_and_WT_M_E18_CORT$cell_type)))
plotDataList$M_MUT_and_WT_M_P30_CORT$Age = c(rep(30, times=length(plotDataList$M_MUT_and_WT_M_P30_CORT$cell_type)))
plotDataList$M_MUT_and_WT_M_P60_CORT$Age = c(rep(60, times=length(plotDataList$M_MUT_and_WT_M_P60_CORT$cell_type)))
plotDataList$M_MUT_and_WT_M_P120_CORT$Age = c(rep(120, times=length(plotDataList$M_MUT_and_WT_M_P120_CORT$cell_type)))
plotDataTogether = do.call(rbind, plotDataList)

plotData = plotDataTogether[-c(1, which(grepl("-activated", plotDataTogether$cell_type) == "TRUE"), which(grepl("-not", plotDataTogether$cell_type) == "TRUE")),]

plotData_male_limma <- plotData %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type) 

male_limma <- ggplot(data=plotData_male_limma, aes(x=Age, y=numDEGs, color=cell_type)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = polychrome_palette)+
  geom_point(size = 3) +
  theme_classic() +
  xlab("Age (days)") +
  ggtitle("limmaVoomCC Postnatal male numDEGs")+
  scale_x_continuous(breaks=c(30, 60, 120))+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'none',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 14, face = 'bold', hjust = 1.0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 14, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/numDEGs_male_postnatal_limmaCC.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
######################
## load female data ##
######################
s.obj.names = c("M_MUT_and_WT_F_P30_CORT", "M_MUT_and_WT_F_P60_CORT", "M_MUT_and_WT_F_P150_CORT")
DEGlists = vector(mode="list", length=3)
plotDataList = vector(mode="list", length=3)

names(DEGlists) = s.obj.names
names(plotDataList) = s.obj.names

for(i in s.obj.names) {
  
  setwd(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/DEG_visualization"))
  
  cell_types = list.dirs(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/DEsingle"))
  cell_types = cell_types[-c(1, which(grepl("QC", cell_types)=="TRUE"), which(grepl("plotData", cell_types)=="TRUE"), which(grepl("interactivePlots", cell_types)=="TRUE"))]
  
  DEGlists[[i]] <- lapply(cell_types, function(cellType) {
    openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
  })
  
  cell_types = gsub(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/DEsingle/"), "", cell_types)
  
  names(DEGlists[[i]]) = cell_types
  
  numDEGs = sapply(cell_types, function(cellType) {
    length(which(DEGlists[[i]][[cellType]][,4] < 0.1))
  })
  
  plotDataList[[i]] = data.frame(cell_type = cell_types, numDEGs = as.numeric(numDEGs))
  
}

### Number of DEGs across time
#female
#plotDataList$M_MUT_and_WT_F_E18_CORT$Age = c(rep(-2, times=length(plotDataList$M_MUT_and_WT_F_E18_CORT$cell_type)))
plotDataList$M_MUT_and_WT_F_P30_CORT$Age = c(rep(30, times=length(plotDataList$M_MUT_and_WT_F_P30_CORT$cell_type)))
plotDataList$M_MUT_and_WT_F_P60_CORT$Age = c(rep(60, times=length(plotDataList$M_MUT_and_WT_F_P60_CORT$cell_type)))
plotDataList$M_MUT_and_WT_F_P150_CORT$Age = c(rep(150, times=length(plotDataList$M_MUT_and_WT_F_P150_CORT$cell_type)))
plotDataTogether = do.call(rbind, plotDataList)

plotData = plotDataTogether[-c(1, which(grepl("-activated", plotDataTogether$cell_type) == "TRUE"), which(grepl("-not", plotDataTogether$cell_type) == "TRUE")),]

plotData_female_desingle <- plotData %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type) 

female_desingle <- ggplot(data=plotData_female_desingle, aes(x=Age, y=numDEGs, color=cell_type)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = polychrome_palette)+
  geom_point(size = 3) +
  theme_classic() +
  xlab("Age (days)") +
  ggtitle("DEsingle Postnatal female numDEGs")+
  scale_x_continuous(breaks=c(30, 60, 150))+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'none',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 14, face = 'bold', hjust = 1.0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 14, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/numDEGs_female_postnatal_DEsingle.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

for(i in s.obj.names) {
  
  setwd(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/DEG_visualization"))
  
  cell_types = list.dirs(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/limmaVoomCC"))
  cell_types = cell_types[-c(1, which(grepl("QC", cell_types)=="TRUE"), which(grepl("plotData", cell_types)=="TRUE"), which(grepl("interactivePlots", cell_types)=="TRUE"))]
  
  DEGlists[[i]] <- lapply(cell_types, function(cellType) {
    openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
  })
  
  cell_types = gsub(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/DEG_two_methods/{i}/limmaVoomCC/"), "", cell_types)
  
  names(DEGlists[[i]]) = cell_types
  
  numDEGs = sapply(cell_types, function(cellType) {
    length(which(DEGlists[[i]][[cellType]][,4] < 0.05))
  })
  
  plotDataList[[i]] = data.frame(cell_type = cell_types, numDEGs = as.numeric(numDEGs))
  
}

### Number of DEGs across time
#female
#plotDataList$M_MUT_and_WT_F_E18_CORT$Age = c(rep(-2, times=length(plotDataList$M_MUT_and_WT_F_E18_CORT$cell_type)))
plotDataList$M_MUT_and_WT_F_P30_CORT$Age = c(rep(30, times=length(plotDataList$M_MUT_and_WT_F_P30_CORT$cell_type)))
plotDataList$M_MUT_and_WT_F_P60_CORT$Age = c(rep(60, times=length(plotDataList$M_MUT_and_WT_F_P60_CORT$cell_type)))
plotDataList$M_MUT_and_WT_F_P150_CORT$Age = c(rep(150, times=length(plotDataList$M_MUT_and_WT_F_P150_CORT$cell_type)))
plotDataTogether = do.call(rbind, plotDataList)

plotData = plotDataTogether[-c(1, which(grepl("-activated", plotDataTogether$cell_type) == "TRUE"), which(grepl("-not", plotDataTogether$cell_type) == "TRUE")),]

plotData_female_limma <- plotData %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type) 

female_limma <- ggplot(data=plotData_female_limma, aes(x=Age, y=numDEGs, color=cell_type)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = polychrome_palette)+
  geom_point(size = 3) +
  theme_classic() +
  xlab("Age (days)") +
  ggtitle("limmaVoomCC Postnatal female numDEGs")+
  scale_x_continuous(breaks=c(30, 60, 150))+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'none',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 14, face = 'bold', hjust = 1.0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 14, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/numDEGs_female_postnatal_limma.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
wrap_plots(female_limma, female_desingle, male_limma, male_desingle) +plot_annotation(tag_levels = 'A')
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/numDEGs_4_panels.pdf",
                device = NULL,
                height = 8.5,
                width = 12)
### Get DEGs that are common to all time points ###
#get rid of activated/not
names(DEGlists[[1]])[c(-3,-4, -6, -7, -9, -10, -12, -13, -15, -16, -19, -20, -22, -23, -25, -26, -28, -29)] 

cell_types.list = lapply(1:3, function(x){
  names(DEGlists[[x]])
})

cell_types.list
cell_types_common = Reduce(intersect, cell_types.list)

intersect.list = lapply(cell_types_common, function(x){
  
  l = lapply(1:3, function(y){
    rownames(DEGlists[[y]][[x]])[which(DEGlists[[y]][[x]]$adj.P.Val<0.1)]
  })
  
  names(l) = names(DEGlists)
  
  intersect = Reduce(intersect, l)
  
  return(intersect)
})

names(intersect.list) = cell_types_common[c(-3,-4, -6, -7, -9, -10, -12, -13, -15, -16, -19, -20, -22, -24, -25, -27, -28)] 

# Cells with persistent differences in AC149090.1 expression
cell_list_0 = c("L2_3_IT", "L2_3_IT-activated", "L2_3_IT-not", "L4-not", "L6-activated")
cell_list_0 <- x
# Cells with persistent differences in Junb expression
cell_list_1 = c("L5-activated", "L6-activated")
cell_list_1 <- x
logFC = sapply(1:3, function(x){
  
  sapply(cell_list_0, function(y){
    DEGlists[[x]][[y]]$logFC[which(rownames(DEGlists[[x]][[y]])=="AC149090.1")]
  })
  
})

logFC.new = logFC %>%
  as_tibble()

logFC.new$cell_type = rownames(logFC)

logFC.new = logFC.new %>%
  gather(key=Age_tmp, value=logFC, V1:V3)

logFC.new$Age = ifelse(logFC.new$Age_tmp=="V1", 30, (ifelse(logFC.new$Age_tmp=="V2", 60, 120)))  

ggplot(logFC.new, aes(x=Age, y=logFC, color=cell_type)) +
  geom_line() +
  theme_classic() +
  xlab("Age (days)") +
  ylab("log(Fold Change) MUT/WT") +
  scale_x_continuous(breaks=c(30, 60, 120)) +
  geom_hline(yintercept=0) +
  ggtitle("AC149090.1")

ggsave("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/LogFC_AC149090.1_Across_Time_Males.pdf")

logFC = sapply(1:3, function(x){
  
  sapply(cell_list_1, function(y){
    DEGlists[[x]][[y]]$logFC[which(rownames(DEGlists[[x]][[y]])=="Junb")]
  })
  
})

logFC.new = logFC %>%
  as_tibble()

logFC.new$cell_type = rownames(logFC)

logFC.new = logFC.new %>%
  gather(key=Age_tmp, value=logFC, V1:V3)

logFC.new$Age = ifelse(logFC.new$Age_tmp=="V1", 30, (ifelse(logFC.new$Age_tmp=="V2", 60, 120)))  

ggplot(logFC.new, aes(x=Age, y=logFC, color=cell_type)) +
  geom_line() +
  theme_classic() +
  xlab("Age (days)") +
  ylab("log(Fold Change) MUT/WT") +
  scale_x_continuous(breaks=c(30, 60, 120)) +
  geom_hline(yintercept=0) +
  ggtitle("Junb")

ggsave("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/LogFC_Junb_Across_Time_Males.pdf")

