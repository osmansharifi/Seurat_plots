#####################
## plot Kegg terms ##
#####################

# By Osman Sharifi 
library(glue)
library(tidyr)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)
library(tidyverse)

###############
## load data ##
###############

master_df <- read.csv("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/master_enricher_GOterms.csv") 

######################################
## set path and order of cell types ##
######################################

pdf_path = "/Users/osman/Desktop/LaSalle_lab/Seurat_figures/"
x = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
master_df<- master_df %>%
  mutate(Cell.Type =  factor(Cell.Type, levels = x)) %>%
  arrange(Cell.Type) 

######################
## Subset master_df ##
######################

#subset E18 males
e18_male <- master_df[which(master_df$Metadata=='M_MUT_and_WT_M_E18_WB' & master_df$Time.Point =="E18"),]
e18_male <- filter(e18_male, Adjusted.P.value <= 0.05)
e18_male <- e18_male %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
e18_male <- e18_male[order(e18_male$Adjusted.P.value),]

#subset E18 females
e18_female <- master_df[which(master_df$Metadata=='M_MUT_and_WT_F_E18_WB' & master_df$Time.Point =="E18"),]
e18_female <- filter(e18_female, Adjusted.P.value <= 0.05)
e18_female <- e18_female %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
e18_female <- e18_female[order(e18_female$Adjusted.P.value),]

#subset P30 females
p30_female <- master_df[which(master_df$Metadata=='M_MUT_and_WT_F_P30_CORT' & master_df$Time.Point =="P30"),]
p30_female <- filter(p30_female, Adjusted.P.value <= 0.05)
p30_female <- p30_female %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p30_female <- p30_female[order(p30_female$Odds.Ratio),]

#subset P30 males
p30_male <- master_df[which(master_df$Metadata=='M_MUT_and_WT_M_P30_CORT' & master_df$Time.Point =="P30"),]
p30_male <- filter(p30_male, Adjusted.P.value <= 0.05)
p30_male <- p30_male %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p30_male <- p30_male[order(p30_male$Odds.Ratio),]

#subset P60 females
p60_female <- master_df[which(master_df$Metadata=='M_MUT_and_WT_F_P60_CORT' & master_df$Time.Point =="P60"),]
p60_female <- filter(p60_female, Adjusted.P.value <= 0.05)
p60_female <- p60_female %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p60_female <- p60_female[order(p60_female$Odds.Ratio),]

#subset P60 males
p60_male <- master_df[which(master_df$Metadata=='M_MUT_and_WT_M_P60_CORT' & master_df$Time.Point =="P60"),]
p60_male <- filter(p60_male, Adjusted.P.value <= 0.05)
p60_male <- p60_male %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p60_male <- p60_male[order(p60_male$Odds.Ratio),]

#subset P150 females
p150_female <- master_df[which(master_df$Metadata=='M_MUT_and_WT_F_P150_CORT' & master_df$Time.Point =="P150"),]
p150_female <- filter(p150_female, Adjusted.P.value <= 0.05)
p150_female <- p150_female %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p150_female <- p150_female[order(p150_female$Odds.Ratio),]

#subset P120 males
p120_male <- master_df[which(master_df$Metadata=='M_MUT_and_WT_M_P120_CORT' & master_df$Time.Point =="P120"),]
p120_male <- filter(p120_male, Adjusted.P.value <= 0.05)
p120_male <- p120_male %>% arrange(Adjusted.P.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p120_male <- p120_male[order(p120_male$Odds.Ratio),]

#########################
## Visualize KEGG Data ##
#########################
#E18 males
E18_male_GO <- ggplot(e18_male,
       aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant E18 Male '
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
ggsave(glue::glue("{pdf_path}E18_Male_KeggTerms_dotplot.pdf"), width = 15,
       height = 12)

#E18 Females
E18_female_GO <- ggplot(e18_female,
                      aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant E18 Female '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}E18_Female_KeggTerms_dotplot.pdf"), width = 15,
       height = 12)

#P30 Females
P30_female_GO <- ggplot(p30_female,
                        aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant P30 Female '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}P30_Female_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)

#P30 males
P30_male_GO <- ggplot(p30_male,
                        aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant P30 male '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}P30_male_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)

#P60 Females
P60_female_GO <- ggplot(p60_female,
                        aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant P60 Female '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}P60_Female_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)

#P60 males
P60_male_GO <- ggplot(p60_male,
                      aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant P60 male '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}P60_male_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)

#P150 Females
P150_female_GO <- ggplot(p150_female,
                        aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant P150 Female '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}P150_Female_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)

#P120 males
P120_male_GO <- ggplot(p120_male,
                      aes(x = Term, y = Cell.Type, size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 KEGG Terms',
    subtitle = 'Significant P120 male '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}P120_male_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)

#########################################################
## Find what females have in common across timepoints ##
#########################################################
p30_col1 <- p30_female$Term
p60_col1 <- p60_female$Term
p150_col1 <- p150_female$Term
length_min <- min(length(p30_col1),length(p60_col1),length(p150_col1))
length_max <- max(length(p30_col1),length(p60_col1),length(p150_col1))
kegg_terms <- data.frame(cbind(p30_col1[1:length_min], p60_col1[1:length_min], p150_col1[1:length_min]))
colnames(kegg_terms)<-c("p30_KEGG","p60_KEGG","p150_KEGG")
common_terms <- data.frame(Reduce(dplyr::intersect, list(kegg_terms$p30_KEGG,kegg_terms$p60_KEGG,kegg_terms$p150_KEGG)))

write.csv(common_terms, glue::glue("{pdf_path}common_terms_only.csv"), row.names = FALSE)
write.csv(p30_female, glue::glue("{pdf_path}p30_female_terms.csv"), row.names = FALSE)
write.csv(p60_female, glue::glue("{pdf_path}p60_female_terms.csv"), row.names = FALSE)
write.csv(p150_female, glue::glue("{pdf_path}p1500_female_terms.csv"), row.names = FALSE)
female_total_kegg <- read.csv(glue::glue("{pdf_path}common_female_terms.csv"))
library(plyr)
female_total_kegg$Time.Point <- revalue(female_total_kegg$Time.Point, c("E18" = 18, "P30" = 30, "P60"= 60, "P120" = 120, "P150" = 150))
female_total_kegg$Time.Point <- as.numeric(female_total_kegg$Time.Point)
female_total_kegg<- female_total_kegg %>%
  mutate(Cell.Type =  factor(Cell.Type, levels = x)) %>%
  arrange(Cell.Type) 

#common female terms plot
female_common <- ggplot(female_total_kegg,
                        aes(x = Term, y = Cell.Type, size = Time.Point, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_size_continuous(breaks = c(30, 60, 150)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Common KEGG Terms',
    subtitle = 'Across time in Females '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}common_female_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)

##########################################################
## Find what the males have in common across timepoints ##
##########################################################
p30_col1_male <- p30_male$Term
p60_col1_male <- p60_male$Term
p120_col1_male <- p120_male$Term
length_min <- min(length(p30_col1_male),length(p60_col1_male),length(p120_col1_male))
length_max <- max(length(p30_col1_male),length(p60_col1_male),length(p120_col1_male))
kegg_terms <- data.frame(cbind(p30_col1_male[1:length_min], p60_col1_male[1:length_min], p120_col1_male[1:length_min]))
colnames(kegg_terms)<-c("p30_KEGG","p60_KEGG","p120_KEGG")
common_terms <- data.frame(Reduce(dplyr::intersect, list(kegg_terms$p30_KEGG,kegg_terms$p60_KEGG,kegg_terms$p120_KEGG)))

write.csv(common_terms, glue::glue("{pdf_path}common_terms_male.csv"), row.names = FALSE)
write.csv(p30_male, glue::glue("{pdf_path}p30_male_terms.csv"), row.names = FALSE)
write.csv(p60_male, glue::glue("{pdf_path}p60_male_terms.csv"), row.names = FALSE)
write.csv(p120_male, glue::glue("{pdf_path}p120_male_terms.csv"), row.names = FALSE)
female_total_kegg <- read.csv(glue::glue("{pdf_path}common_female_terms.csv"))
library(plyr)
female_total_kegg$Time.Point <- revalue(female_total_kegg$Time.Point, c("E18" = 18, "P30" = 30, "P60"= 60, "P120" = 120, "P150" = 150))
female_total_kegg$Time.Point <- as.numeric(female_total_kegg$Time.Point)
female_total_kegg<- female_total_kegg %>%
  mutate(Cell.Type =  factor(Cell.Type, levels = x)) %>%
  arrange(Cell.Type) 

#common female terms plot
female_common <- ggplot(female_total_kegg,
                        aes(x = Term, y = Cell.Type, size = Time.Point, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_size_continuous(breaks = c(30, 60, 150)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Common KEGG Terms',
    subtitle = 'Across time in Females '
  )  +   
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
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
ggsave(glue::glue("{pdf_path}common_female_KEGGTerms_dotplot.pdf"), width = 15,
       height = 12)


