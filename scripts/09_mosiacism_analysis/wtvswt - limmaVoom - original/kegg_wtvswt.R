# create a common kegg term plot between different timepoints of mouse WT vs WT

library(dplyr)
library(ggplot2)
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
# Load data
kegg_terms <- readxl::read_xlsx(glue::glue("{base_path}broad_group_analysis/GABA_kegg_all.xlsx",col_names = TRUE))

kegg_terms <- filter(kegg_terms, P.value <= 0.05)
kegg_terms$Time_Point = kegg_terms$Time_point
# Get the top 5 terms for each cell type
top_terms <- kegg_terms %>%
  group_by(Time_Point) %>%
  top_n(10, P.value) 

top_terms$Time_Point <- factor(top_terms$Time_Point, levels = c("P30", "P60", "P150"))
ggplot(top_terms, aes(x = Term, y = Time_Point, color = P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_gradientn(name = "P-value", colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu")) +
  scale_size_continuous(range = c(3, 8)) +
  scale_y_discrete(name = "Time point") +
  scale_x_discrete(name = "Terms") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'Top10 GABAergic KEGG terms') +
  theme(legend.position = "bottom")+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 12, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 8, face = 'plain', vjust = 0.5, colour = "black"),
    axis.title = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.x = element_text(size = 18, face = 'bold', colour = "black"),
    axis.title.y = element_text(size = 18, face = 'bold', colour = "black"),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"), # Text size
    title = element_text(size = 18, face = "bold")) +
  coord_flip()
ggsave(glue::glue("{base_path}mut_from_mut_vs_wt_from_wt/GABA_top10_kegg.pdf"), 
       device = NULL,
       height = 8.5,
       width = 7)

###########################################
## Venn diagram of the overlapping genes ##
###########################################
library(VennDiagram)
library(grDevices)
# Extract significant DEGs
sig_kegg_P30 <- rownames(deg_results$P30)[deg_results$P30$adj.P.Val < 0.05]
sig_kegg_P60 <- rownames(deg_results$P60)[deg_results$P60$adj.P.Val < 0.05]
sig_kegg_P150 <- rownames(deg_results$P150)[deg_results$P150$adj.P.Val < 0.05]
intersection_all2 <- intersect(sig_genes_P150,sig_genes_P30)
intersection_all3 <- intersect(intersection_all2,sig_genes_P60)
# Extract KEGG terms
p30_rows <- subset(kegg_terms, Time_Point == "P30")
# Filter rows for P60
p60_rows <- subset(kegg_terms, Time_Point == "P60")
# Filter rows for P150
p150_rows <- subset(kegg_terms, Time_Point == "P150")
intersection_all2 <- intersect(p30_rows$Term,p60_rows$Term)
intersection_all3 <- intersect(intersection_all2,p150_rows$Term)
# Create a Venn diagram
pdf(file=glue("{base_path}mut_from_mut_vs_wt_from_wt/venn_kegg_GABA.pdf"))
temp<- venn.diagram(
  x = list(
    P30 = p30_rows$Term,
    P60 = p60_rows$Term,
    P150 = p150_rows$Term
  ),
  category.names = c("P30", "P60", "P150"),
  main = 'GABA KEGG Terms of MUT cells from MUT females and WT cells from WT females ',
  filename = NULL,
  col = c("#E6B8BFFF","#CC7A88FF","#990F26FF"),
  fill = c("#E6B8BFFF","#CC7A88FF","#990F26FF"),
  cat.cex = 1.2,
  cat.fontface = "bold",
  euler.d = TRUE,
  disable.logging = TRUE,
  hyper.test = TRUE
)
grid.draw(temp)
dev.off()
# Write the merged dataframe to a CSV file
write.csv(merged, file = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_male_mouse_kegg_overlap.csv", row.names = FALSE)