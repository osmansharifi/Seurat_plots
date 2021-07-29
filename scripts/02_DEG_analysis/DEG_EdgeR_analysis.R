library(Seurat)
library(edgeR)
library(SingleCellExperiment)
library(scran)
library(glue)

# By Viktoria Haghani

################################################################################
## Variables

# Paths

# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_E18_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P30_with_labels_proportions.rda"
data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P60_with_labels_proportions.rda"
# data_file <- "~/GitHub/snRNA-seq-pipeline/raw_data/rett_P120_with_labels_proportions3.rda"

DEG_data_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_M_P60_CORT"

################################################################################
## Data preparation

# Load in data
load(data_file)
experiment.aggregate
Idents(experiment.aggregate) <- 'celltype.call'

# Prepare data
# Rename "Non-neuronal" as "Non_neuronal" for variable name usage
experiment.aggregate <- RenameIdents(object = experiment.aggregate, "Non-neuronal" = "Non_neuronal")
# We want to get rid of the G2M and S phase cells, so subset to keep only G1 cells
experiment.aggregate <- subset(x = experiment.aggregate, subset = cell.cycle == "G1")
# Set mitochondrial threshold to 0.5%
experiment.aggregate <- subset(x = experiment.aggregate, subset = percent.mito <= "0.5")

################################################################################
## EdgeR Analysis

for (cell_type in cell_types){
  # Subset Seurat object for the cluster
  edgeR_obj <- subset(experiment.aggregate, subset = celltype.call == cell_type)
  # Generate a count matrix 
  counts <- as.matrix(edgeR_obj@assays$RNA@counts)
  counts <- counts[Matrix::rowSums(counts >= 1) >= 1, ]
  # Subset the meta data for filtered gene/cells
  metadata <- edgeR_obj@meta.data
  metadata <- metadata[,c("nCount_RNA", "orig.ident")]
  metadata <- metadata[colnames(counts),]
  # Make single cell experiment
  sce <- SingleCellExperiment(assays = counts, 
                              colData = metadata)
  # Convert to edgeR object
  dge <- convertTo(sce, type="edgeR", assay.type = 1)
  meta_dge <- dge$samples
  meta_dge <- meta_dge[,c("lib.size","norm.factors")]
  meta_dge <- cbind(meta_dge, metadata)
  meta_dge$group <- factor(meta_dge$orig.ident)
  levels(meta_dge$group)
  meta_dge$group <- relevel(meta_dge$group, "MUT_M_P60_CORT1", "MUT_M_P60_CORT2", "WT_M_P60_CORT1", "WT_M_P60_CORT2")
  dge$samples <- meta_dge
  # Model fit
  dge <- calcNormFactors(dge)
  design <- model.matrix(~0+group, data=dge$samples)
  head(design)
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  # Differential expression testing
  my.contrasts <- makeContrasts(MUT_vs_WT = c(groupWT_M_P60_CORT1+groupWT_M_P60_CORT2) - c(groupMUT_M_P60_CORT1+groupMUT_M_P60_CORT2), levels = design)
  qlf.contrast <- glmQLFTest(fit, contrast=my.contrasts)
  # Use Benjamini-Hochberg correction for p-values
  qlf.contrast.all.genes <- topTags(qlf.contrast, n = 1000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
  # Note that topTags() outputs adjusted p-values in the FDR column
  cell_EdgeR_DEG <- subset(x = qlf.contrast.all.genes$table, subset = FDR < 0.05)
  # Write data to CSV so analysis does not need to be rerun when working with data
  write.csv(cell_EdgeR_DEG, file = glue(DEG_data_dir, cell_type, "_", metadata_info, "_EdgeR_DEG.csv"))
}
