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
  theme_minimal()
