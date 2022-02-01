library(ggplot2)
library(dplyr)
library(viridis)
# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths
#csv_path <- "~/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/master_go_data_males.csv"
#figure_path <- "~/GitHub/snRNA-seq-pipeline/figures/go_analysis/"
csv_path <- "~/Documents/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/master_go_data_top3.csv"
figure_path <- "~/Documents/GitHub/snRNA-seq-pipeline/figures/go_analysis/"
# Names
plot_title <- "GO Analysis"
plot_subtitle <- "All Cell Types, Time Points, and Ontologies for Male Mice"

################################################################################

# Read in GO Data from master CSV file
go_data <- read.csv(file = csv_path)

#subset and sort GO terms
E18_GO <- go_data[which(go_data$Time.Point=='E18'),]
E18_GO <- E18_GO[order(E18_GO$Fisher),]
p.adjust <- round(p.adjust(E18_GO$Fisher,method="BH"),digits = 4)
E18_GO=cbind(E18_GO,p.adjust)
E18_GO <- E18_GO[order(E18_GO$p.adjust),]
E18_GO <- filter(E18_GO, p.adjust < 0.05)
E18_GO$GeneRatio <- E18_GO$Significant/E18_GO$Annotated

P30_GO <- go_data[which(go_data$Time.Point=='P30'),]
P30_GO <- P30_GO[order(P30_GO$Fisher),]
p.adjust <- round(p.adjust(P30_GO$Fisher,method="BH"),digits = 4)
P30_GO=cbind(P30_GO,p.adjust)
P30_GO <- P30_GO[order(P30_GO$p.adjust),]
P30_GO <- filter(P30_GO, p.adjust < 0.05)
P30_GO$GeneRatio <- P30_GO$Significant/P30_GO$Annotated

P60_GO <- go_data[which(go_data$Time.Point=='P60'),]
P60_GO <- P60_GO[order(P60_GO$Fisher),]
p.adjust <- round(p.adjust(P60_GO$Fisher,method="BH"),digits = 4)
P60_GO=cbind(P60_GO,p.adjust)
P60_GO <- P60_GO[order(P60_GO$p.adjust),]
P60_GO <- filter(P60_GO, p.adjust < 0.05)
P60_GO$GeneRatio <- P60_GO$Significant/P60_GO$Annotated
P60_GO <- P60_GO[order(P60_GO$GeneRatio),]

P120_GO <- go_data[which(go_data$Time.Point=='P120'),]
P120_GO <- P120_GO[order(P120_GO$Fisher),]
p.adjust <- round(p.adjust(P120_GO$Fisher,method="BH"),digits = 4)
P120_GO=cbind(P120_GO,p.adjust)
P120_GO <- P120_GO[order(P120_GO$p.adjust),]
P120_GO <- filter(P120_GO, p.adjust < 0.05)
P120_GO$GeneRatio <- P120_GO$Significant/P120_GO$Annotated


# Visualize GO Data
#########################################################################################################
gg1 <- ggplot(E18_GO,
              aes(x = Term, y = Cell.Type, size = GeneRatio, fill = p.adjust))  +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'FDR corrected TopGO Terms',
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
ggsave("E18_Male_GOTerms_dotplot.pdf",gg1, width = 15,
       height = 12)
#P30
gg2 <- ggplot(P30_GO,
              aes(x = Term, y = Cell.Type, size = GeneRatio, fill = p.adjust))  +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'FDR corrected TopGO Terms',
    subtitle = 'Significant P30 Male'
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
ggsave("P30_Male_GOTerms_dotplot.pdf",gg2, width = 15,
       height = 12)

#P60
gg3 <- ggplot(P60_GO,
              aes(x = Term, y = Cell.Type, size = GeneRatio, fill = p.adjust))  +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'FDR corrected TopGO Terms',
    subtitle = 'Significant P60 Male'
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
ggsave("P60_Male_GOTerms_dotplot.pdf",gg3, width = 15,
       height = 12)

#P120
gg4 <- ggplot(P120_GO,
              aes(x = Term, y = Cell.Type, size = GeneRatio, fill = p.adjust))  +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'FDR corrected TopGO Terms',
    subtitle = 'Significant P120 Male'
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
ggsave("P120_Male_GOTerms_dotplot.pdf",gg4, width = 15,
       height = 12)
