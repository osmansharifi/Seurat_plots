library(Seurat)
library(cowplot)
library(ggplot2)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths
 data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_E18_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P30_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P60_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P120_with_labels_proportions.rda"
# data_file <- "/Users/osman/Desktop/LaSalle_lab/Scripts/P30_script/P30_Male_Cortex/rett_P30_with_labels_proportions.rda"

 data_vis_dir <- "~/GitHub/snRNA-seq-pipeline/figures/data_structure_visualization/M_MUT_and_WT_M_E18_WB"
# data_vis_dir <- "~/GitHub/snRNA-seq-pipeline/figures/data_structure_visualization/M_MUT_and_WT_M_P30_CORT"
# data_vis_dir <- "~/GitHub/snRNA-seq-pipeline/figures/data_structure_visualization/M_MUT_and_WT_M_P60_CORT"
# data_vis_dir <- "~/GitHub/snRNA-seq-pipeline/figures/data_structure_visualization/M_MUT_and_WT_M_P120_CORT"
# data_vis_dir <- "/Users/osman/Desktop/LaSalle_lab/Scripts/P30_script/P30_Male_Cortex/reanalyzed"


# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "Mice, Male, E18, WB"

################################################################################
## Visualization and data preparation

# Load in data
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'

# Values represent cell numbers for each cell type
before_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident)

# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, 'Non-neuronal' = 'Non_neuronal')

## Subset cells in G1 and visualize UMAP

# Visualize UMAP
pcs.use <- 10
experiment.aggregate <- RunUMAP(experiment.aggregate, dims = 1:pcs.use)

# This will show us total cells including those in G2M and S phase
cell_type_grouping_including_G2M_and_S <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "cell.cycle") +
  ggtitle("Cell Type Grouping Including G2M and S Phase Cells", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_type_grouping_including_G2M_and_S.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

# We want to get rid of the G2M and S phase cells, so subset to keep only G1 cells
experiment.aggregate <- subset(x = experiment.aggregate, subset = cell.cycle == "G1")

# Generating a UMAP plot to validate that G2M and S phase cells were removed
cell_type_grouping_for_subsetted_data_G1_only <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "cell.cycle") +
  ggtitle("Cell Type Grouping for Subsetted Data (G1 Only)", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_type_grouping_for_subsetted_data_G1_only.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

# Generate UMAP plot with cell types
cell_types_after_G1_subsetting <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "celltype.call", label = TRUE) +
  ggtitle("Cell Types after G1 Subsetting", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_types_after_G1_subsetting.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

# Visualize cell type per sample
cell_types_by_orig_ident_before_subset <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Cell Types by Experimental and Control Groups after G1 Subsetting", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_types_by_orig_ident_before_subset.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

## Subset to remove mitochondrial genes

# See percent mitochondrial genes before subsetting to threshold
percent_mt_before_subset <- VlnPlot(experiment.aggregate, features = "percent.mito") +
  ggtitle("% mt gene expression before subsetting", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14)) +
  ylab("Expression (%)")
ggsave("percent_mt_before_subset.pdf", device = "pdf", path = data_vis_dir, width = 15, height = 9)

# Set threshold to 0.5%
experiment.aggregate <- subset(x = experiment.aggregate, subset = percent.mito <= "0.5")

# Validate the mitochondrial genes are removed
percent_mt_after_subset <- VlnPlot(experiment.aggregate, features = "percent.mito") +
  ggtitle("% mt gene expression after setting 0.5%  threshold", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14)) +
  ylab("Expression (%)")
ggsave("percent_mt_after_subset.pdf", device = "pdf", path = data_vis_dir, width = 15, height = 9)

# Values represent cell numbers for each cell type
after_subset_cell_counts <- table(Idents(experiment.aggregate), experiment.aggregate$orig.ident)

## See changes in cell number after subsetting

# Call table to compare cells before G2M, S, and mt exclusion
before_subset_cell_counts
# Generate a data frame from the table
before_subset_cell_counts_df <- data.frame(before_subset_cell_counts)
# Order values so bars appear in descending value order
before_subset_cell_counts_df$Var1 <- reorder(before_subset_cell_counts_df$Var1,-before_subset_cell_counts_df$Freq)
sample_name <- before_subset_cell_counts_df$Var2
# Create grouped bar plot
cell_counts_before_subset <- ggplot(before_subset_cell_counts_df, aes(fill=sample_name, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Cell Counts Per Condition Before Subsetting", subtitle = metadata_info) +
  xlab("Cell Types") +
  ylab("Number of Cells") +
  ylim(0, 12500) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("cell_counts_before_subset.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

# Call table to compare cells left after  G2M, S, and mt exclusion
after_subset_cell_counts 
# Generate a data frame from the table
after_subset_cell_counts_df <- data.frame(after_subset_cell_counts)
# Order values so bars appear in descending value order
after_subset_cell_counts_df$Var1 <- reorder(after_subset_cell_counts_df$Var1,-after_subset_cell_counts_df$Freq)
sample_name <- after_subset_cell_counts_df$Var2
# Create grouped bar plot
cell_counts_after_subset <- ggplot(after_subset_cell_counts_df, aes(fill=sample_name, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Cell Counts Per Condition After Subsetting", subtitle = metadata_info) +
  xlab("Cell Types") +
  ylab("Number of Cells") +
  ylim(0, 500) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("cell_counts_after_subset.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

## Visualize clusters after G2M, S, and mt exclusion

# Generate UMAP plot with cell types after G1 and mitochondrial subsetting
cell_types_after_subsetting_G1_mt_subset <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "celltype.call", label = TRUE) +
  ggtitle("Cell Types After G1 & mt Subsetting", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_types_after_subsetting_G1_mt_subset.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

# Visualize UMAP plot grouped by experimental adn control groups to see how cell types match each group
cell_type_by_orig_ident_after_subset <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Cell Types per Group After G1 & mt Subsetting", subtitle = metadata_info)  +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("cell_type_by_orig_ident_after_subset.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)

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
combined_replicates_after_subsetting <- DimPlot(experiment.aggregate, reduction = "umap", group.by = "new.ident") +
  ggtitle("Cell Types per Group After Combining MUT and WT", subtitle = metadata_info) +
  theme(legend.text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 14))
ggsave("combined_replicates_after_subsetting.pdf", device = "pdf", path = data_vis_dir, width = 10, height = 9)