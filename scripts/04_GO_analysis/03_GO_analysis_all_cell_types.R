library(ggplot2)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths
#csv_path <- "~/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/master_go_data_males.csv"
#figure_path <- "~/GitHub/snRNA-seq-pipeline/figures/go_analysis/"
csv_path <- "~/Documents/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/master_go_data_males.csv"
figure_path <- "~/Documents/GitHub/snRNA-seq-pipeline/figures/go_analysis/"
# Names
plot_title <- "GO Analysis"
plot_subtitle <- "All Cell Types, Time Points, and Ontologies for Male Mice"

################################################################################

# Read in GO Data from master CSV file
go_data <- read.csv(file = csv_path)

# Visualize GO Data

ggplot(go_data, aes(x = Metadata, y = Term, size = -log10(go_data[9]))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  xlab('') +
  ylab('Enrichment score') +
  labs(title = plot_title, subtitle = plot_subtitle) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1.5, 3)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    # Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  coord_flip()

#ggplot2::ggsave(glue('{figure_path}master_go_data_Fisher.pdf'),
#                device = NULL,
#                height = 8.5,
#                width = 12)