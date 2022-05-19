##Author: Osman Sharifi
library(Seurat) 
library(dplyr)
library(ggplot2)
library(plyr)

# Custom color palette
polychrome_palette  <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")
# Order that celltypes should appear in 
x = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
# Load data
load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all.cortex.combined.RData")
meta.data <- all.cortex.combined[[]]

# create random classifications for the sake of this example
meta.data$cell_type <- sample(meta.data$celltype.call, nrow(meta.data), replace = TRUE)
meta.data <- meta.data %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type) 

counts <- group_by(meta.data, cell_type, Age, Sex, Condition) %>% dplyr::summarise(count = n())
counts$Age <- revalue(counts$Age, c("E18" = 18, "P30" = 30, "P60"= 60, "P120" = 120, "P150" = 150))
counts$Age <- as.numeric(counts$Age)
#counts$Age <- counts %>% dplyr::arrange(Age)
counts_postnatal_wt_male <- filter(counts, Sex =="Male", Age!=18, Condition=="WT")
counts_postnatal_mut_male <- filter(counts, Sex =="Male", Age!=18, Condition=="MUTANT")
counts_postnatal_wt_female <- filter(counts, Sex =="Female", Age!=18, Condition=="WT")
counts_postnatal_mut_female <- filter(counts, Sex =="Female", Age!=18, Condition=="MUTANT")
counts_E18 <- filter(counts, Age!= 30, Age!=60, Age!=120, Age!=150)

write.csv(counts_postnatal_mut_male, glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/counts_postnatal_mut_male.csv"), row.names = FALSE)
write.csv(counts_postnatal_wt_male, glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/counts_postnatal_wt_male.csv"), row.names = FALSE)
write.csv(counts_postnatal_mut_female, glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/counts_postnatal_mut_female.csv"), row.names = FALSE)
write.csv(counts_postnatal_wt_female, glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/counts_postnatal_wt_female.csv"), row.names = FALSE)
write.csv(counts_E18, glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/counts_E18.csv"), row.names = FALSE)

##plot the line graphs
wt_m_plot <- ggplot(counts_postnatal_wt_male, aes(Age, count, color = cell_type), group = 3) +
  geom_line()+
  geom_point()+
  theme_classic() +
  xlab("Age (days)") +
  ylab("numCells") +
  ggtitle("numCells_WT_males")+
  scale_color_manual(values = polychrome_palette)+
  scale_x_continuous(breaks=c(30, 60, 120))

mut_m_plot <- ggplot(counts_postnatal_mut_male, aes(Age, count, color = cell_type), group = 3) +
  geom_line(size = 2)+
  geom_point()+
  theme_classic() +
  xlab("Age (days)") +
  ylab("propCells") +
  ggtitle("proportion of cells in mutant males")+
  scale_color_manual(values = polychrome_palette)+
  scale_x_continuous(breaks=c(30, 60, 120))

wt_f_plot <- ggplot(counts_postnatal_wt_female, aes(Age, count, color = cell_type), group = 3) +
  geom_line()+
  geom_point()+
  theme_classic() +
  xlab("Age (days)") +
  ylab("numCells") +
  ggtitle("numCells_WT_females")+
  scale_color_manual(values = polychrome_palette)+
  scale_x_continuous(breaks=c(30, 60, 150))

mut_f_plot <- ggplot(counts_postnatal_mut_female, aes(Age, count, color = cell_type), group = 3) +
  geom_line()+
  geom_point()+
  theme_classic() +
  xlab("Age (days)") +
  ylab("numCells") +
  ggtitle("numCells_MUT_females")+
  scale_color_manual(values = polychrome_palette)+
  scale_x_continuous(breaks=c(30, 60, 150))

wrap_plots(A. = wt_m_plot, B.=mut_m_plot, C. = wt_f_plot, D.=mut_f_plot)
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/numCells.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

ggplot(data=counts_E18, aes(x=cell_type, y=count, fill = Condition)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal()+
  scale_fill_brewer(palette="Set1")+
  theme(axis.text.x = element_text(angle = 45))+
  xlab("Cell_type") +
  ylab("numCells") +
  ggtitle("numCells_E18") +
  coord_polar()

ggplot(counts_E18, aes(x=Condition, y=count, fill=cell_type)) +
  geom_bar(stat="identity", width=1, color="white")# +
 # coord_polar("y", start=0)# +
  #theme_void() + #
 # theme(legend.position="none") +
  #geom_text(aes(label = cell_type), color = "white", size=6) +
  #scale_fill_brewer(palette="Set1")
  
####################################
## proportions of cells over time ##
####################################
prop_postnatal_mut_male <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/prop_postnatal_mut_male.csv", header = TRUE)
prop_postnatal_mut_female <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/prop_postnatal_mut_female.csv", header = TRUE)
prop_postnatal_wt_male <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/prop_postnatal_wt_male.csv", header = TRUE)
prop_postnatal_wt_female <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/prop_postnatal_wt_female.csv", header = TRUE)
prop_E18 <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/cell_data/prop_E18.csv", header = TRUE)
prop_postnatal_mut_male <- prop_postnatal_mut_male %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type)
prop_postnatal_mut_female <- prop_postnatal_mut_female %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type)
prop_postnatal_wt_male <- prop_postnatal_wt_male %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type)
prop_postnatal_wt_female <- prop_postnatal_wt_female %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type)
prop_E18 <- prop_E18 %>%
  mutate(cell_type =  factor(cell_type, levels = x)) %>%
  arrange(cell_type)

######################
## plot proportions ##
######################
wt_m_plot_prop <- ggplot(prop_postnatal_wt_male, aes(Age, count, color = cell_type), group = 3) +
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  theme_classic() +
  xlab("Age (days)") +
  ylab("propCells") +
  ggtitle("propCells_WT_males")+
  scale_color_manual(values = polychrome_palette)+
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

mut_m_plot_prop <- ggplot(prop_postnatal_mut_male, aes(Age, count, color = cell_type), group = 3) +
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  theme_classic() +
  xlab("Age (days)") +
  ylab("propCells") +
  ggtitle("propCells_MUT_males")+
  scale_color_manual(values = polychrome_palette)+
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

wt_f_plot_prop <- ggplot(prop_postnatal_wt_female, aes(Age, count, color = cell_type), group = 3) +
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  theme_classic() +
  xlab("Age (days)") +
  ylab("propCells") +
  ggtitle("propCells_WT_females")+
  scale_color_manual(values = polychrome_palette)+
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

mut_f_plot_prop <- ggplot(prop_postnatal_mut_female, aes(Age, count, color = cell_type), group = 3) +
  geom_line(size = 1.5)+
  geom_point(size = 3)+
  theme_classic() +
  xlab("Age (days)") +
  ylab("propCells") +
  ggtitle("propCells_MUT_females")+
  scale_color_manual(values = polychrome_palette)+
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

wrap_plots(A. = wt_m_plot_prop, B.=mut_m_plot_prop, C. = wt_f_plot_prop, D.=mut_f_plot_prop)+plot_annotation(tag_levels = 'A')
ggplot2::ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/propCells.pdf",
                device = NULL,
                height = 8.5,
                width = 12)

