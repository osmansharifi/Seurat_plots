# library
library(ggplot2)
library(data.table)
library(magrittr)
library(viridis)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)

## extract meta data
md <- human_rettcort@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
test <- md[, .N, by = c("Samples", "celltype")]

## with additional casting after the counting
md[, .N, by = c("Samples", "predicted.subclass")] %>% dcast(., Samples ~ predicted.subclass, value.var = "N")

#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

# Stacked + percent
ggplot(test, aes(fill=celltype, y=N, x=Samples)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values= polychrome_palette)+
  ylab("Proportions") +
labs(
  title = 'Cell type proportions',
  subtitle = 'Human cortices'
)  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
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
    title = element_text(size = 14, face = "bold")) +
  coord_flip()

# Best way to get cell type proportions
# Order of samples
experiment.merged$orig.ident <- factor(experiment.merged$orig.ident, levels=c("WT_F_E18_WB1", "WT_F_E18_WB2","MUT_F_E18_WB1", "MUT_F_E18_WB2", "WT_M_E18_WB1", "WT_M_E18_WB2", "MUT_M_E18_WB1", "MUT_M_E18_WB2"))

# Order of cell types
experiment.merged$celltype <- factor(experiment.merged$celltype, levels=c("L2_3_IT", "L5", "L5_6_NP", "L6", "Pvalb", "Vip", "Sst", "Sncg", "Lamp5", "Peri", "Endo", "Oligo", "Astro", "Non-neuronal", "Meis2", "OPC", "VLMC"))
#plot the cell proportions
experiment.merged@meta.data %>%
  group_by(orig.ident,celltype) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=orig.ident,y=percent, fill=celltype)) +
  geom_col() +
  ggtitle("Percentage of cell types in the embryonic brain") +
  scale_fill_manual(values = polychrome_palette)+
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

