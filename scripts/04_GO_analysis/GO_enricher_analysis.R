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

master_df_female <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/04_GO_analysis/female_GO.csv") 
master_df_male <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/04_GO_analysis/male_GO.csv") 
######################################
## set path and order of cell types ##
######################################

pdf_path = "/Users/osman/Desktop/LaSalle_lab/Seurat_figures/"
x = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
master_df_female<- master_df_female %>%
  mutate(Cell.Type =  factor(Cell.Type, levels = x)) %>%
  arrange(Cell.Type) 
master_df_male<- master_df_male %>%
  mutate(Cell.Type =  factor(Cell.Type, levels = x)) %>%
  arrange(Cell.Type) 

######################
## Subset master_df ##
######################

#subset E18 males
e18_male <- master_df_male[which(master_df_male$Metadata=='M_MUT_and_WT_M_E18_WB' & master_df_male$Time.Point =="E18"),]
e18_male <- filter(e18_male, X.log10.p.value >= 3.3)
e18_male <- e18_male %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
e18_male <- e18_male[order(e18_male$X.log10.p.value),]

#subset E18 females
e18_female <- master_df_female[which(master_df_female$Metadata=='M_MUT_and_WT_F_E18_WB' & master_df_female$Time.Point =="E18"),]
e18_female <- filter(e18_female, X.log10.p.value >= 3.3)
e18_female <- e18_female %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
e18_female <- e18_female[order(e18_female$X.log10.p.value),]

#subset P30 females
p30_female <- master_df[which(master_df$Metadata=='M_MUT_and_WT_F_P30_CORT' & master_df$Time.Point =="P30"),]
p30_female <- filter(p30_female, X.log10.p.value <= 0.05)
p30_female <- p30_female %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p30_female <- p30_female[order(p30_female$Odds.Ratio),]

p30_female_limma <- master_df_female[which(master_df_female$Metadata=='M_MUT_and_WT_F_P30_CORT' & master_df_female$Time.Point =="P30" & master_df_female$deg_method == "LimmaVoomCC"),]
p30_female_limma <- filter(p30_female_limma, X.log10.p.value >= 3.3)
p30_female_limma <- p30_female_limma %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p30_female_limma <- p30_female_limma[order(p30_female_limma$X.log10.p.value),]

#subset P30 males
p30_male_limma <- master_df_male[which(master_df_male$Metadata=='M_MUT_and_WT_M_P30_CORT' & master_df_male$Time.Point =="P30" & master_df_male$deg_method == "LimmaVoomCC"),]
p30_male_limma <- filter(p30_male_limma, X.log10.p.value >= 3.3)
p30_male_limma <- p30_male_limma %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p30_male_limma <- p30_male_limma[order(p30_male_limma$X.log10.p.value),]

#subset P60 females
p60_female_limma <- master_df_female[which(master_df_female$Metadata=='M_MUT_and_WT_F_P60_CORT' & master_df_female$Time.Point =="P60" & master_df_female$deg_method == "LimmaVoomCC"),]
p60_female_limma <- filter(p60_female_limma, X.log10.p.value >= 3.3)
p60_female_limma <- p60_female_limma %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p60_female_limma <- p60_female_limma[order(p60_female_limma$X.log10.p.value),]

#subset P60 males
p60_male_limma <- master_df_male[which(master_df_male$Metadata=='M_MUT_and_WT_M_P60_CORT' & master_df_male$Time.Point =="P60" & master_df_male$deg_method == "LimmaVoomCC"),]
p60_male_limma <- filter(p60_male_limma, X.log10.p.value >= 3.3)
p60_male_limma <- p60_male_limma %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p60_male_limma <- p60_male_limma[order(p60_male_limma$X.log10.p.value),]

#subset P150 females
p150_female_limma <- master_df_female[which(master_df_female$Metadata=='M_MUT_and_WT_F_P150_CORT' & master_df_female$Time.Point =="P150" & master_df_female$deg_method == "LimmaVoomCC"),]
p150_female_limma <- filter(p150_female_limma, X.log10.p.value >= 3.3)
p150_female_limma <- p150_female_limma %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p150_female_limma <- p150_female_limma[order(p150_female_limma$X.log10.p.value),]

#subset P120 males
p120_male_limma <- master_df_male[which(master_df_male$Metadata=='M_MUT_and_WT_M_P120_CORT' & master_df_male$Time.Point =="P120" & master_df_male$deg_method == "LimmaVoomCC"),]
p120_male_limma <- filter(p120_male_limma, X.log10.p.value >= 3.3)
p120_male_limma <- p120_male_limma %>% arrange(X.log10.p.value) %>%
  group_by(Cell.Type) %>% top_n(10)
p120_male_limma <- p120_male_limma[order(p120_male_limma$X.log10.p.value),]

#########################
## Visualize KEGG Data ##
#########################
#E18 males
E18_male_GO <- ggplot(e18_male,
       aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis(option = "plasma", direction = -1) + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}E18_Male_GOTerms.pdf"), width = 15,
       height = 12)

#E18 Females
E18_female_GO <- ggplot(e18_female,
                      aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}E18_Female_GOTerms.pdf"), width = 15,
       height = 12)

#P30 Females
P30_female_GO <- ggplot(p30_female_limma,
                        aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}P30_Female_GOTerms.pdf"), width = 15,
       height = 12)

#P30 males
P30_male_GO <- ggplot(p30_male_limma,
                        aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}P30_male_GOTerms.pdf"), width = 15,
       height = 12)

#P60 Females
P60_female_GO <- ggplot(p60_female_limma,
                        aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}P60_Female_GOTerms.pdf"), width = 15,
       height = 12)

#P60 males
P60_male_GO <- ggplot(p60_male_limma,
                      aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}P60_male_GOTerms.pdf"), width = 15,
       height = 12)

#P150 Females
p150_female_GO <- ggplot(p150_female_limma,
                        aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}p150_female_GOTerms.pdf"), width = 15,
       height = 12)

#P120 males
P120_male_GO <- ggplot(p120_male_limma,
                      aes(x = Term, y = Cell.Type, size = X.log10.p.value, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Top 10 GO Terms',
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
ggsave(glue::glue("{pdf_path}P120_male_GOTerms.pdf"), width = 15,
       height = 12)

#########################################################
## Find what females have in common across timepoints ##
#########################################################
p30_col1 <- p30_female_limma$Term
p60_col1 <- p60_female_limma$Term
p150_col1 <- p150_female_limma$Term
length_min <- min(length(p30_col1),length(p60_col1),length(p150_col1))
length_max <- max(length(p30_col1),length(p60_col1),length(p150_col1))
go_terms <- data.frame(cbind(p30_col1[1:length_min], p60_col1[1:length_min], p150_col1[1:length_min]))
colnames(go_terms)<-c("p30_GO","p60_GO","p150_GO")
common_terms <- data.frame(Reduce(dplyr::intersect, list(go_terms$p30_GO,go_terms$p60_GO,go_terms$p150_GO)))

write.csv(common_terms, glue::glue("{pdf_path}common_female_GOterms_only.csv"), row.names = FALSE)
write.csv(p30_female_limma, glue::glue("{pdf_path}p30_female_GOterms.csv"), row.names = FALSE)
write.csv(p60_female_limma, glue::glue("{pdf_path}p60_female_GOterms.csv"), row.names = FALSE)
write.csv(p150_female_limma, glue::glue("{pdf_path}p150_female_GOterms.csv"), row.names = FALSE)
female_total_go <- read.csv(glue::glue("{pdf_path}common_female_GOterms.csv"))
library(plyr)
female_total_go$Time.Point <- revalue(female_total_go$Time.Point, c("E18" = 18, "P30" = 30, "P60"= 60, "P120" = 120, "P150" = 150))
female_total_go$Time.Point <- as.numeric(female_total_go$Time.Point)
female_total_go<- female_total_go %>%
  mutate(Cell.Type =  factor(Cell.Type, levels = x)) %>%
  arrange(Cell.Type) 
female_go_limma <- female_total_go[which(female_total_go$deg_method == "LimmaVoomCC"),]

#common female terms plot
female_common <- ggplot(female_total_go,
                        aes(x = Term, y = Cell.Type, size = Time.Point, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_size_continuous(breaks = c(30, 60, 150)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Common GO Terms',
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
ggsave(glue::glue("{pdf_path}common_female_GOTerms.pdf"), width = 15,
       height = 12)

##########################################################
## Find what the males have in common across timepoints ##
##########################################################
p30_col1_male <- p30_male_limma$Term
p60_col1_male <- p60_male_limma$Term
p120_col1_male <- p120_male_limma$Term
length_min <- min(length(p30_col1_male),length(p60_col1_male),length(p120_col1_male))
length_max <- max(length(p30_col1_male),length(p60_col1_male),length(p120_col1_male))
go_terms <- data.frame(cbind(p30_col1_male[1:length_min], p60_col1_male[1:length_min], p120_col1_male[1:length_min]))
colnames(go_terms)<-c("p30_GO","p60_GO","p120_GO")
common_terms <- data.frame(Reduce(dplyr::intersect, list(go_terms$p30_GO,go_terms$p60_GO,go_terms$p120_GO)))

write.csv(common_terms, glue::glue("{pdf_path}common_GOterms_male_only.csv"), row.names = FALSE)
write.csv(p30_male_limma, glue::glue("{pdf_path}p30_male_GOterms.csv"), row.names = FALSE)
write.csv(p60_male_limma, glue::glue("{pdf_path}p60_male_GOterms.csv"), row.names = FALSE)
write.csv(p120_male_limma, glue::glue("{pdf_path}p120_male_GOterms.csv"), row.names = FALSE)
male_total_go <- read.csv(glue::glue("{pdf_path}common_GOterms_male.csv"))
library(plyr)
male_total_go$Time.Point <- revalue(male_total_go$Time.Point, c("E18" = 18, "P30" = 30, "P60"= 60, "P120" = 120, "P150" = 150))
male_total_go$Time.Point <- as.numeric(male_total_go$Time.Point)
male_total_go<- male_total_go %>%
  mutate(Cell.Type =  factor(Cell.Type, levels = x)) %>%
  arrange(Cell.Type) 
male_go_limma <- male_total_go[which(male_total_go$deg_method == "LimmaVoomCC"),]

#common female terms plot
male_common <- ggplot(male_total_go,
                        aes(x = Term, y = Cell.Type, size = Time.Point, fill = X.log10.p.value)) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_size_continuous(breaks = c(30, 60, 120)) +
  scale_fill_viridis() + 
  xlab('') + ylab('Cell Type') +
  labs(
    title = 'Common GO Terms',
    subtitle = 'Across time in males '
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
ggsave(glue::glue("{pdf_path}common_male_GOTerms.pdf"), width = 15,
       height = 12)

wrap_plots(A. = male_common, B.=female_common)+plot_annotation(tag_levels = 'A')
