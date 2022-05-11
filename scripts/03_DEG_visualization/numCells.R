##Author: Osman Sharifi
library(Seurat) 
library(dplyr)
library(ggplot2)
library(plyr)

#Custom color palette
polychrome_palette  <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")
# Order that celltypes should appear in 
x = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
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
  geom_line()+
  geom_point()+
  theme_classic() +
  xlab("Age (days)") +
  ylab("numCells") +
  ggtitle("numCells_MUT_males")+
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
