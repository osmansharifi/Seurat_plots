# create a common kegg term plot between human and mouse

library(dplyr)
library(ggplot2)
# Load human
human <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_total_kegg.csv", header = TRUE)
human$Sex <- "females"
human$Specie <- "Human"
human <- filter(human, P.value <= 0.05)
# Load mouse
mouse <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/mouse_total_kegg.csv", header = TRUE)
female_mouse <- filter(mouse, Sex =="females")
female_mouse$Specie <- "Mouse"
female_mouse <- filter(female_mouse, P.value <= 0.05)
# Merge the human and mouse data frames by their common columns
merged <- merge(human, female_mouse, by = c("Term", "Cell_Type"))

# Get the top 5 terms for each cell type
top_terms <- merged %>%
  group_by(Cell_Type) %>%
  top_n(5, P.value.x) 

# Print the top 5 terms for each cell type
print(top_terms)

ggplot(top_terms, aes(x = Term, y = Cell_Type, color = P.value.x)) +
  geom_point(size = 8) +
  scale_color_gradientn(name = "P-value", colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu")) +
  scale_size_continuous(range = c(3, 8)) +
  scale_y_discrete(name = "Cell Type") +
  scale_x_discrete(name = "Terms") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'Top5 overlapping KEGG terms between Human and Female Mouse') +
  theme(legend.position = "bottom")+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 18, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 18, face = 'bold', vjust = 0.5, colour = "black"),
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
ggsave(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_female_mouse_kegg_overlap.tiff"), width = 15,
       height = 12)
# Write the merged dataframe to a CSV file
write.csv(merged, file = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_female_mouse_kegg_overlap.csv", row.names = FALSE)

# Merge the human and mouse data frames by their common columns
male_mouse <- filter(mouse, Sex =="males")
male_mouse <- filter(male_mouse, P.value <= 0.05)
merged <- merge(human, male_mouse, by = c("Term", "Cell_Type"))

# Get the top 5 terms for each cell type
top_terms <- merged %>%
  group_by(Cell_Type) %>%
  top_n(5, P.value.x) 

ggplot(top_terms, aes(x = Term, y = Cell_Type, color = P.value.x)) +
  geom_point(size = 8) +
  scale_color_gradientn(name = "P-value", colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu")) +
  scale_size_continuous(range = c(3, 8)) +
  scale_y_discrete(name = "Cell Type") +
  scale_x_discrete(name = "Terms") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'Top5 overlapping KEGG terms between Human and Male Mouse') +
  theme(legend.position = "bottom")+
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 90, size = 18, face = 'bold', hjust = 1.0, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 18, face = 'bold', vjust = 0.5, colour = "black"),
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
ggsave(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_male_mouse_kegg_overlap.tiff"), 
       width = 16,
       height = 12)
# Write the merged dataframe to a CSV file
write.csv(merged, file = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_male_mouse_kegg_overlap.csv", row.names = FALSE)
