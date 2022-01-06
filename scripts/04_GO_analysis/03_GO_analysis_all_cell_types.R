library(ggplot2)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths
csv_path <- "~/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/master_go_data_males.csv"
figure_path <- "~/GitHub/snRNA-seq-pipeline/figures/go_analysis/enrichment_scores/"

# Names
plot_title <- "GO Analysis"
plot_subtitle <- "All Cell Types, Time Points, and Ontologies for Male Mice"

################################################################################

# Read in GO Data from master CSV file
go_data <- read.csv(file = csv_path)

# Visualize GO Data

ggplot(data = go_data, aes(x = Metadata, y = Term, 
                        color = Fisher, size = Significant)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("Cell Types for Each Time Point") + 
  labs(title = plot_title, subtitle = plot_subtitle)