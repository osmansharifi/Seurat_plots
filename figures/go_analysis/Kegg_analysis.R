#####################
## plot Kegg terms ##
#####################

# By Osman Sharifi 
library(glue)
library(tidyr)
library(ggplot2)
library(viridis)
library(dplyr)

###############
## load data ##
###############
master_df <- read.csv("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/master_enricher_GOterms.csv") 

pdf_path = "/Users/osman/Desktop/LaSalle_lab/Seurat_figures/"

######################
## Subset master_df ##
######################

#subset and sort GO terms
E18_GO <- male_go[which(male_go$Metadata=='M_MUT_and_WT_M_E18_WB' & male_go$Time.Point =="E18"),]
E18_GO <- E18_GO[order(E18_GO$`-log10(p.value)`),]
E18_GO <- filter(E18_GO, `-log10(p.value)` >= 1.3)
E18_GO <- E18_GO %>% arrange(desc(`-log10(p.value)`)) %>%
  group_by(Gene.Ontology) %>% top_n(40)

#subset and sort GO terms female E18
E18_GO_Female <- female_go[which(female_go$Metadata=='M_MUT_and_WT_F_E18_WB' & female_go$Time.Point =="E18"),]
E18_GO_Female <- E18_GO_Female[order(E18_GO_Female$`-log10(p.value)`),]
E18_GO_Female <- filter(E18_GO_Female, `-log10(p.value)` >= 1.3)
E18_GO_Female <- E18_GO_Female %>% arrange(desc(`-log10(p.value)`)) %>%
  group_by(Gene.Ontology) %>% top_n(40)

P30_GO_Female <- female_go[which(female_go$Metadata=='M_MUT_and_WT_F_P30_CORT' & female_go$Time.Point =="P30"),]
P30_GO_Female <- P30_GO_Female[order(P30_GO_Female$`-log10(p.value)`),]
P30_GO_Female <- filter(P30_GO_Female, `-log10(p.value)` >= 1.3)
P30_GO_Female <- P30_GO_Female %>% arrange(desc(`-log10(p.value)`)) %>%
  group_by(Gene.Ontology) %>% top_n(50)

P30_GO_male <- male_go[which(male_go$Metadata=='M_MUT_and_WT_M_P30_CORT' & male_go$Time.Point =="P30"),]
P30_GO_male <- P30_GO_male[order(P30_GO_male$`-log10(p.value)`),]
P30_GO_male <- filter(P30_GO_male, `-log10(p.value)` >= 1.3)
P30_GO_male <- P30_GO_male %>% arrange(desc(`-log10(p.value)`)) %>%
  group_by(Gene.Ontology) %>% top_n(40)

P60_GO <- go_data[which(go_data$Time.Point=='P60'),]
P60_GO <- P60_GO[order(P60_GO$Fisher),]
p.adjust <- round(p.adjust(P60_GO$Fisher,method="BH"),digits = 4)
P60_GO=cbind(P60_GO,p.adjust)
P60_GO <- P60_GO[order(P60_GO$p.adjust),]
P60_GO <- filter(P60_GO, p.adjust < 0.05)
P60_GO$GeneRatio <- P60_GO$Significant/P60_GO$Annotated
P60_GO <- P60_GO[order(P60_GO$GeneRatio),]

#subset and sort GO terms
P60_GO_Female <- go_data[which(go_data$Time.Point=='P60'),]
P60_GO_Female <- P60_GO_Female[which(P60_GO_Female$Sex=='F'),]
P60_GO_Female <- P60_GO_Female[order(P60_GO_Female$Fisher),]
p.adjust <- round(p.adjust(P60_GO_Female$Fisher,method="BH"),digits = 4)
P60_GO_Female=cbind(P60_GO_Female,p.adjust)
P60_GO_Female <- P60_GO_Female[order(P60_GO_Female$p.adjust),]
P60_GO_Female <- filter(P60_GO_Female, p.adjust < 0.05)
P60_GO_Female$GeneRatio <- P60_GO_Female$Significant/P60_GO_Female$Annotated

P120_GO <- go_data[which(go_data$Time.Point=='P120'),]
P120_GO <- P120_GO[order(P120_GO$Fisher),]
p.adjust <- round(p.adjust(P120_GO$Fisher,method="BH"),digits = 4)
P120_GO=cbind(P120_GO,p.adjust)
P120_GO <- P120_GO[order(P120_GO$p.adjust),]
P120_GO <- filter(P120_GO, p.adjust < 0.05)
P120_GO$GeneRatio <- P120_GO$Significant/P120_GO$Annotated

#subset and sort GO terms
P150_GO_Female <- go_data[which(go_data$Time.Point=='P150'),]
P150_GO_Female <- P150_GO_Female[which(P150_GO_Female$Sex=='F'),]
P150_GO_Female <- P150_GO_Female[order(P150_GO_Female$Fisher),]
p.adjust <- round(p.adjust(P150_GO_Female$Fisher,method="BH"),digits = 4)
P150_GO_Female=cbind(P150_GO_Female,p.adjust)
P150_GO_Female <- P150_GO_Female[order(P150_GO_Female$p.adjust),]
P150_GO_Female <- filter(P150_GO_Female, p.adjust < 0.05)
P150_GO_Female$GeneRatio <- P150_GO_Female$Significant/P150_GO_Female$Annotated


# Visualize GO Data
#########################################################################################################
E18_male_GO <- ggplot(E18_GO,
       aes(x = Term, y = Cell.Type, size = `-log10(p.value)`, fill = `-log10(p.value)`))  +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top GO Terms',
    subtitle = 'Significant E18 Male '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 12, face = 'bold', hjust = 1.0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  coord_flip()
ggsave(glue::glue("{pdf_path}E18_Male_GOTerms_dotplot.pdf"), width = 15,
       height = 12)

#E18 Females

E18_female_GO <- ggplot(E18_GO_Female,
                      aes(x = Term, y = Cell.Type, size = `-log10(p.value)`, fill = `-log10(p.value)`))  +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top GO Terms',
    subtitle = 'Significant E18 Female '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 12, face = 'bold', hjust = 1.0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  coord_flip()
ggsave(glue::glue("{pdf_path}E18_Female_GOTerms_dotplot.pdf"), width = 15,
       height = 12)

#P30 females

P30_female_GO <- ggplot(P30_GO_Female,
                        aes(x = Term, y = Cell.Type, size = `-log10(p.value)`, fill = `-log10(p.value)`))  +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top GO Terms',
    subtitle = 'Significant P30 Female '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 12, face = 'bold', hjust = 1.0, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  coord_flip()
ggsave(glue::glue("{pdf_path}P30_Female_GOTerms_dotplot.pdf"), width = 15,
       height = 12)

master_ericher <- read.csv("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/master_enricher_GOterms.csv")  
