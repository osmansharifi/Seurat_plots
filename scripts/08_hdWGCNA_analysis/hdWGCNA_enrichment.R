# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(1234)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- load('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/mouse_cortex_WGCNA.RData')

# enrichr databases to test
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023', 'KEGG_2019_Mouse')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  adult_postnatal,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "KEGG_2019_Mouse", # this has to be one of the lists we used above!!!
  n_terms=3 # number of terms for each module
)
ggplot2::ggsave(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/kegg_module.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)
# Assuming enrich_df is your DataFrame
# If not, replace it with the actual DataFrame variable name
# For example: enrich_df <- data.frame(Column1 = c(1, 2, 3), Column2 = c('A', 'B', 'C'))

# Specify the directory path
directory_path <- "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/"

# Create the CSV file path
csv_file_path <- paste0(directory_path, "enrich_df.csv")

# Save the DataFrame to CSV
write.csv(enrich_df, file = csv_file_path, row.names = FALSE)

cat("DataFrame has been saved to", csv_file_path, "\n")

# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- 'celltype.call'
markers <- Seurat::FindAllMarkers(
  adult_postnatal,
  only.pos = TRUE,
  logfc.threshold=1)

DefaultAssay(seurat_obj) <- 'RNA'
markers <- FindAllMarkers(object = seurat_obj) %>%
  Add_Pct_Diff() %>%
  filter(pct_diff > 0.6)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
wrap_plots(plot_list, ncol=3)

OverlapDotPlot <- function(
    overlap_df, plot_var = 'odds_ratio',
    logscale=TRUE,
    neglog=FALSE,
    plot_significance=TRUE,
    ...
){
  
  label <- plot_var
  if(logscale){
    overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
    label <- paste0('log(', plot_var, ')')
  }
  if(neglog){
    overlap_df[[plot_var]] <- -1 * overlap_df[[plot_var]]
    label <- paste0('-', label)
  }
  
  p <- overlap_df %>% ggplot(aes(x=module, y=group)) +
    geom_point(aes(
      size=get(plot_var)),
      #alpha=get(plot_var)),
      color=overlap_df$color
    ) +
    RotatedAxis() +
    ylab('') + xlab('') + labs(size=label) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  # plot significance level?
  if(plot_significance){
    p <- p + geom_text(aes(label=Significance))
  }
  
  p
}
# plot odds ratio of the overlap as a dot plot
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio',logscale = TRUE) +
  ggtitle('Overlap of modules & cell-type markers') +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 10) +
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
    
    # Legend
    legend.key = element_blank(),  # removes the border
    legend.key.size = unit(1, "cm"),  # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"),  # Text size
    title = element_text(size = 18, face = "bold"),
  )
ggplot2::ggsave(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/module_marker_overlap.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)
csv_file_path <- paste0(directory_path, "marker_module_overlap.csv")
write.csv(overlap_df, file = csv_file_path, row.names = FALSE)
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 6,
  text_digits = 6,
  text_color = 'black',
  high_color = '#B2182B',
  mid_color = '#EEEEEE',
  low_color = '#2166AC',
  plot_max = 0.2,
  combine=FALSE
)[1] #this specifies the celltype
ggplot2::ggsave(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/all_cells_module_trait.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

# Extract the relevant information from meta.data
plot_data <- seurat_obj@meta.data[, c("Sex", "turquoise", "Age", "Genotype")]
plot_data$Age <- factor(plot_data$Age, levels = c("P30", "P60", "P120", "P150"))
# Create a violin plot with different colors for male and female
p <- ggplot(plot_data, aes(x = Age, y = turquoise, fill = Sex)) +
  geom_violin() +
  labs(title = "Violin Plot of Time point vs Turquoise",
       x = "Time point",
       y = "Module Turquoise Eigennode") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'Module Trait correlation') +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 10) +
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
    
    # Legend
    legend.key = element_blank(),  # removes the border
    legend.key.size = unit(1, "cm"),  # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"),  # Text size
    title = element_text(size = 18, face = "bold"),
  )
# Add p-value annotation
p + stat_compare_means(comparisons = list(c("P30", "P60"), c("P30", "P120"), c("P30", "P150"), c("P60", "P120"), c("P60", "P150"), c("P120", "P150")), method = "t.test", label = "p.format")
ggplot2::ggsave(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/Timepoint_eigennode.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)


# Assuming plot_data is your data frame
# You can use the t.test function to compare turquoise values between the two sexes for each age category
# Create an empty data frame to store the results
results_df <- data.frame(Age_Category = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Get unique age categories
unique_age_categories <- unique(plot_data$Age)

# Loop through each age category
for (age_category in unique_age_categories) {
  
  # Subset data for the current age category
  subset_data <- subset(plot_data, Age == age_category)
  
  # Check the age category
  if (age_category == "P120") {
    # For P120, compare it to P150
    t_test_result <- t.test(turquoise ~ 1, data = subset_data, subset = (Age %in% c("P120", "P150")))
  } else if (age_category == "P150") {
    # Skip the comparison for P150
    next
  } else {
    # For other age categories, compare 'Male' to 'Female'
    t_test_result <- t.test(turquoise ~ Sex, data = subset_data, subset = (Sex %in% c("Male", "Female")))
  }
  
  # Store the result in the data frame
  results_df <- rbind(results_df, data.frame(Age_Category = age_category, p_value = t_test_result$p.value))
}
csv_file_path <- paste0(directory_path, "trait_module_statistics.csv")
write.csv(results_df, file = csv_file_path, row.names = FALSE)

# Create a violin plot with different colors for genotype
p <- ggplot(plot_data, aes(x = Genotype, y = turquoise, fill = Sex)) +
  geom_violin() +
  labs(title = "Violin Plot of Genotype vs Turquoise",
       x = "Genotype",
       y = "Module Turquoise Eigennode") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'Module Trait correlation') +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 10) +
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
    
    # Legend
    legend.key = element_blank(),  # removes the border
    legend.key.size = unit(1, "cm"),  # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"),  # Text size
    title = element_text(size = 18, face = "bold"),
  )

# Add p-value annotation
p + stat_compare_means(comparisons = list(c("WT", "MUT")), method = "t.test", label = ('p.format'))
ggplot2::ggsave(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/Genotype_eigennode.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

mut <- filtered_plot_data <- plot_data %>%
  filter(Genotype == "MUT")

# Create a violin plot with different colors for genotype
p <- ggplot(mut, aes(x = Sex, y = turquoise, fill = Sex)) +
  geom_violin() +
  labs(title = "Violin Plot of Genotype vs Turquoise",
       x = "Genotype",
       y = "Module Turquoise Eigennode") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  labs(title = 'Module Trait correlation') +
  theme(legend.position = "bottom") +
  theme_bw(base_size = 10) +
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
    
    # Legend
    legend.key = element_blank(),  # removes the border
    legend.key.size = unit(1, "cm"),  # Sets overall area/size of the legend
    legend.text = element_text(size = 18, face = "bold"),  # Text size
    title = element_text(size = 18, face = "bold"),
  )

# Add p-value annotation
p + stat_compare_means(comparisons = list(c("Male", "Female")), method = "t.test", label = ('p.format'))
ggplot2::ggsave(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/Genotype_eigennode.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)