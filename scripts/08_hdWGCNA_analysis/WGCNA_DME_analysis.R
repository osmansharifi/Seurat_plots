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
set.seed(12345)
setwd('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis')

# re-load the Zhou et al snRNA-seq dataset processed with hdWGCNA
seurat_obj <- load('mouse_cortex_postnatal.RData')