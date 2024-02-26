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

# plot odds ratio of the overlap as a dot plot
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
ggplot2::ggsave(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/module_marker_overlap.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

test <- subset(seurat_obj, subset = seurat_obj@meta.data[['predicted_cell_type']] == '0:CD8 T cell')

PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'black',
  high_color = '#B2182B',
  mid_color = '#EEEEEE',
  low_color = '#2166AC',
  plot_max = 0.2,
  combine=FALSE
)[1] #this specifies the celltype

