# This program will perform differential module eigengene (DME) analysis on single cell data
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrepel)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(1234)
setwd('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis')

# re-load the Zhou et al snRNA-seq dataset processed with hdWGCNA
adult_postnatal <- load('mouse_cortex_postnatal.RData')

# DME analysis comparing WT to RETT 
WT <- adult_postnatal@meta.data %>% subset(Genotype == 'WT') %>% colnames()[,1]
MUT <- adult_postnatal@meta.data %>% subset(Genotype == 'MUT') %>% rownames

DMEs <- FindDMEs(
  adult_postnatal,
  barcodes1 = WT,
  barcodes2 = MUT,
  test.use='wilcox',
  wgcna_name='Genotype'
)

head(DMEs)

Idents(adult_postnatal) <- "Genotype"
test <- CalculateBarcodeInflections(adult_postnatal)
BarcodeInflectionsPlot(test) + scale_colour_manual(name = "orig.ident", values = polychrome_palette)

test2 <- rownames(test)
test2 <- gsub("*-", "",test2)

# create a new dataframe with the split row names
test_split <- data.frame(do.call(rbind, strsplit(row.names(test), "-")))

# name the columns
colnames(test_split) <- c("Barcodes", "Sample_name")

# split the Barcodes column by _ and keep only the info to the right of the _
test_split$Barcodes <- sub(".*_", "", test_split$Barcodes)

# merge the test_split and adult_postnatal@meta.data dataframes
merged <- merge(test_split, adult_postnatal@meta.data, by.x="Sample_name", by.y="orig.ident", all.x=TRUE)

# replace the existing Barcodes column in adult_postnatal@meta.data with the new one from merged
adult_postnatal@meta.data$Barcodes <- merged$Barcodes
