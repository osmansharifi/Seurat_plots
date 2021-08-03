library(Seurat)
library(cowplot)
library(ggplot2)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths
data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_E18_with_labels_proportions.rda"

data_vis_dir <- "~/GitHub/snRNA-seq-pipeline/figures/data_structure_visualization/M_MUT_and_WT_M_E18_WB"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "Mice, Male, E18, Whole Brain"

################################################################################
## Visualization 

# No G1 and mt subsetting for embryonic because most cells are cycling

# Load in data
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'

# Values represent cell numbers for each cell type
cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident)

# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')

## Visualize UMAP

# Visualize UMAP
pcs.use <- 10
experiment.aggregate <- RunUMAP(experiment.aggregate, dims = 1:pcs.use)

# Generate UMAP plot with cell cycle phases
cell_cycle_phase_grouping <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "cell.cycle") +
  ggtitle("Cell Cycle Phase Grouping", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_cycle_phase_grouping.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

# Generate UMAP plot with cell types
cell_types <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "celltype.call", label = TRUE) +
  ggtitle("Cell Types", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_types.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

# Visualize cell type per sample by identity
cell_types_by_orig_ident <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Cell Types by Experimental and Control Groups", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_types_by_orig_ident.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

percent_mt <- VlnPlot(experiment.aggregate, features = "percent.mito") +
  ggtitle("% mt gene expression", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14)) +
  ylab("Expression (%)")
ggsave("percent_mt.pdf", device = "pdf", path = data_vis_dir, width = 15, height = 9)

## Look at cell_counts

# Generate a data frame from the table
cell_counts_df <- data.frame(cell_counts)
# Order values so bars appear in descending value order
cell_counts_df$Var1 <- reorder(cell_counts_df$Var1,-cell_counts_df$Freq)
sample_name <- cell_counts_df$Var2
# Create grouped bar plot
cell_counts <- ggplot(cell_counts_df, aes(fill=sample_name, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Cell Counts Per Condition", subtitle = metadata_info) +
  xlab("Cell Types") +
  ylab("Number of Cells") +
  ylim(0, 12500) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("cell_counts.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

## Reorganize Seurat object identities for DEG analysis

# Create only MUT and WT groups
experiment.aggregate@meta.data$new.ident <- plyr::mapvalues(
  x = experiment.aggregate@meta.data$orig.ident, 
  from = c("MUT_M_E18_WB1", "MUT_M_E18_WB2", "WT_M_E18_WB1", "WT_M_E18_WB2"), 
  to = c("MUT_M_E18_WB", "MUT_M_E18_WB", "WT_M_E18_WB", "WT_M_E18_WB")
)

## See counts for orig.ident and validate that they're combined correctly for new.idents
old_ident_counts <- table(experiment.aggregate@meta.data$orig.ident)
old_ident_counts
new_ident_counts <- table(experiment.aggregate@meta.data$new.ident)
new_ident_counts

# See grouping after to validate proper separation and data structure
combined_replicates <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "new.ident") +
  ggtitle("Cell Types per Group After Combining MUT and WT", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("combined_replicates.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

