library(Seurat)
library(dplyr)
library(scCustomize)
library(patchwork)

# load parsed table
counts_table <- read.csv('/Users/osman/Desktop/LaSalle_lab/Rett_Data/E18/cell_parsing_outputs/cell_parsing_df_bodyadded.csv')
adult_male <- filter(counts_table, Sex == 'M' & Timepoint != 'E18')
table(adult_male$Sample, adult_male$Timepoint)
adult_male$X <- NULL
adult_male$Timepoint <- NULL
adult_male$Sex <- NULL
adult_male$Condition <- NULL
adult_male$Sample <- NULL

#Load data
load("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/p120.RData")
Idents(p120) <- "celltype.call"
male.all = p120
levels(male.all) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
DefaultAssay(p120) <- 'RNA'
head(male.all)

# Extract the barcodes from rownames using regular expressions
new_rownames <- sub(".*_(.*?)-.*", "\\1", rownames(male.all@meta.data))

# Add a unique identifier to each rowname
unique_rownames <- make.unique(new_rownames, sep = "_")

# Update the rownames in the dataset
rownames(male.all@meta.data) <- unique_rownames

# putting counts back into seurat
barcode_counts <- adult_male[, c("Barcode", "WT_Mecp2", "MUT_Mecp2")]
rownames(barcode_counts) <- barcode_counts$Barcode
barcode_counts$Barcode <- NULL

# Assuming your_dataframe is the name of your dataframe

# Add a unique identifier to each barcode using an index
unique_barcodes <- make.unique(barcode_counts$Barcode, sep = "_")
index <- seq_along(unique_barcodes)

# Create a new dataframe with the unique barcodes and counts
barcode_counts <- data.frame(Barcode = unique_barcodes, WT_Mecp2 = barcode_counts$WT_Mecp2, MUT_Mecp2 = barcode_counts$MUT_Mecp2)

# Set the rownames using the combination of barcode and index
rownames(barcode_counts) <- paste0(barcode_counts$Barcode, "_", index)
barcode_counts$Barcode <- NULL
# Convert the barcode counts dataframe to a matrix
counts_matrix <- as.matrix(barcode_counts)

# Add the counts matrix as features to the Seurat object
male.all[["WT-Mecp2"]] <- CreateAssayObject(counts_matrix)
male.all[["MUT-Mecp2"]] <- CreateAssayObject(counts_matrix)


library(stringr)

# Extract the barcodes from Seurat object's meta.data
barcodes <- str_extract(rownames(p120@meta.data), "(?<=_).*(?=-)")
# Check if every row in the Barcode column is unique
is_unique <- !any(duplicated(adult_male$Barcode))
# Print the result
print(is_unique)
# Create a new column 'Mecp2_status' in Seurat object's meta.data
p120@meta.data$Mecp2_status <- "None"

# Find matching barcodes and transfer corresponding info
matching_rows <- barcodes %in% adult_male$Barcode
p120@meta.data$Mecp2_status[matching_rows] <- adult_male$Condition[match(barcodes[matching_rows], adult_male$Barcode)]

DimPlot_scCustom(seurat_object = p120, pt.size = 0.6, split.by = 'Condition', group.by = 'Mecp2_status', order = TRUE, label.box = TRUE, ggplot_default_colors = TRUE)

# Subset on a value in the object meta data
mut_test <- subset(x = p120, subset = Mecp2_status == "MUT")
DimPlot_scCustom(seurat_object = p120, pt.size = 0.6, split.by = 'Condition', group.by = 'Mecp2_status', order = TRUE, label.box = TRUE, ggplot_default_colors = TRUE)
cl <- c('gray', 'red')
DimPlot_scCustom(seurat_object = p120, pt.size = 0.8, split.by = 'Condition', group.by = 'Mecp2_allele', order = TRUE, label.box = TRUE, colors_use = c('gray', 'green3'))
ggplot2::ggsave("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/p120_Mecp2allele_green_UMAP.tiff",
                device = NULL,
                height = 10,
                width = 12)
#Visualizations
FeaturePlot(male.all, features = c('WT-Mecp2', 'MUT-Mecp2'), raster = FALSE,pt.size = 0.6, max.cutoff = 10, split.by = 'Condition', order = TRUE)
VlnPlot_scCustom(seurat_object = male.all, features = c('Oaz1', 'Hnrnpa2b1', 'Rpl8', 'Calm1', 'AC149090.1', 'Actg1', 'Pdlim7', 'Snrnp70', 'Plk2', 'Pea15a', 'Map1b', 'Atp6v0c'), split.by = 'Condition', pt.size = 0.9, ggplot_default_colors = TRUE, group.by = 'celltype.call', raster = FALSE, add.noise = FALSE)

VlnPlot_scCustom(seurat_object = male.all, features = c('Oaz1', 'Atp5g3', 'Rpl8', 'AC149090.1'), split.by = 'Condition', pt.size = 0.9, group.by = 'celltype.call', raster = FALSE, ggplot_default_colors = TRUE)
Idents(male.all)
L2_3_IT <- subset(x = male.all, subset = celltype.call == "L2_3_IT")
Idents(L2_3_IT) <- 'Condition'
L2_3_IT.de.markers <- FindMarkers(L2_3_IT, ident.1 = "MUTANT", ident.2 = "WT", logfc.threshold = 0.4)        
VlnPlot_scCustom(seurat_object = L2_3_IT, features = rownames(L2_3_IT.de.markers), split.by = 'Condition', pt.size = 0.9, raster = FALSE, ggplot_default_colors = TRUE, num_columns = 5)      
