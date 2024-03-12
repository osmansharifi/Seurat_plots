library(tidyr)
library(dplyr)
library(ggplot2)

# Read the CSV file with read_csv from readr
data <- read_csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/Mecp2_allele_counts.csv")

# Create a new column 'Gender' based on the row identifier
data <- data %>%
  mutate(Gender = ifelse(grepl("female", ...1, ignore.case = TRUE), "Female", "Male"))

# Reshape the data to long format
data_long <- gather(data, key = "Category", value = "Count", -...1, -Gender)

# Create a stacked bar graph
ggplot(data_long, aes(x = ...1, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Gender, scales = "free_x", space = "free_x", switch = "x") +
  labs(title = "Mecp2 allele count",
       x = "Genotype",
       y = "Cell count") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
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
    title = element_text(size = 18, face = "bold")) 
ggsave(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/Mecp2_allele_counts.pdf"), width = 15,
       height = 12)
