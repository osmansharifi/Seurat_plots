# Required
install.packages("Seurat")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("glue")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("MAST")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")









# Original (I move up the ones I know I need)
install.packages("limma") 
install.packages("S4Vectors")
install.packages("stats4")
install.packages("BiocGenerics")
install.packages("parallel")
install.packages("patchwork")
install.packages("tidyverse")
install.packages("Matrix.utils")
install.packages("dplyr")
install.packages("Matrix")
install.packages("pheatmap")
install.packages("purrr")
install.packages("reshape2")
install.packages("tibble")
install.packages("png")
install.packages("RColorBrewer")
 
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("edgeR")

 
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
 
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("scater")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("apeglm")


# Library load
library(Seurat)
library(ggplot2)
library(edgeR)
library(limma) # edgeR dependency
library(DESeq2)
library(S4Vectors) # DESeq2 dependency
library(stats4) # DESeq2 dependency
library(BiocGenerics) # DESeq2 dependency
library(parallel) # DESeq2 dependency
library(patchwork)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix.utils)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(tibble)
library(pheatmap)
library(apeglm)
library(png)
library(RColorBrewer)
library(MAST)
