library(Seurat)
library(MAST)
library(limma)
library(edgeR)
library(DESeq2)
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(ggVennDiagram)
library(edgeR)
library(udunits2)
library(scran)
library(Rcpp)

# Paths
EdgeR_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/M_MUT_and_WT_M_P60_CORT/"
DESeq2_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/M_MUT_and_WT_M_P60_CORT/"
Limma_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_P60_CORT/"

venn_dir <- "~/GitHub/snRNA-seq-pipeline/figures/venn_diagrams/M_MUT_and_WT_M_P60_CORT/"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_M_P60_CORT"
subtitle_info <- "Mice, Male, P60, Cortex"

################################################################################
# Venn Diagram for Differentially Expressed Genes Per Analysis
# By Viktoria Haghani

for (cell_type in cell_types){
  # Read in data from Limma analysis
  assign(paste0(cell_type, "_Limma_DEG"), read.csv(file = glue(Limma_DEG_dir, cell_type, "_", metadata_info, "_Limma_DEG.csv")))
  # Read in data from EdgeR analysis
  assign(paste0(cell_type, "_EdgeR_DEG"), read.csv(file = glue(EdgeR_DEG_dir, cell_type, "_", metadata_info, "_EdgeR_DEG.csv")))
  # Read in data from DESeq2 analysis
  assign(paste0(cell_type, "_DESeq2_DEG"), read.csv(file = glue(DESeq2_DEG_dir, cell_type, "_", metadata_info, "_DESeq2_DEG.csv")))
}

# List of genes differentially expressed per cluster for Limma
L2_3_IT_Limma_gene_list <- L2_3_IT_Limma_DEG$X
L6_Limma_gene_list <- L6_Limma_DEG$X
Sst_Limma_gene_list <- Sst_Limma_DEG$X
L5_Limma_gene_list <- L5_Limma_DEG$X
L4_Limma_gene_list <- L4_Limma_DEG$X
Pvalb_Limma_gene_list <- Pvalb_Limma_DEG$X
Sncg_Limma_gene_list <- Sncg_Limma_DEG$X
Non_neuronal_Limma_gene_list <- Non_neuronal_Limma_DEG$X
Oligo_Limma_gene_list <- Oligo_Limma_DEG$X
Vip_Limma_gene_list <- Vip_Limma_DEG$X
Lamp5_Limma_gene_list <- Lamp5_Limma_DEG$X
Astro_Limma_gene_list <- Astro_Limma_DEG$X
Peri_Limma_gene_list <- Peri_Limma_DEG$X
Endo_Limma_gene_list <- Endo_Limma_DEG$X
unique_Limma_genes <- unique(c(L2_3_IT_Limma_DEG$X,
                     L6_Limma_DEG$X,
                     Sst_Limma_DEG$X,
                     L5_Limma_DEG$X, 
                     L4_Limma_DEG$X,
                     Pvalb_Limma_DEG$X,
                     Sncg_Limma_DEG$X,
                     Non_neuronal_Limma_DEG$X,
                     Oligo_Limma_DEG$X,
                     Vip_Limma_DEG$X,
                     Lamp5_Limma_DEG$X,
                     Astro_Limma_DEG$X,
                     Peri_Limma_DEG$X,
                     Endo_Limma_DEG$X))

# List of genes differentially expressed per cluster for DESeq2
L2_3_IT_DESeq2_gene_list <- L2_3_IT_DESeq2_DEG$X
L6_DESeq2_gene_list <- L6_DESeq2_DEG$X
Sst_DESeq2_gene_list <- Sst_DESeq2_DEG$X
L5_DESeq2_gene_list <- L5_DESeq2_DEG$X
L4_DESeq2_gene_list <- L4_DESeq2_DEG$X
Pvalb_DESeq2_gene_list <- Pvalb_DESeq2_DEG$X
Sncg_DESeq2_gene_list <- Sncg_DESeq2_DEG$X
Non_neuronal_DESeq2_gene_list <- Non_neuronal_DESeq2_DEG$X
Oligo_DESeq2_gene_list <- Oligo_DESeq2_DEG$X
Vip_DESeq2_gene_list <- Vip_DESeq2_DEG$X
Lamp5_DESeq2_gene_list <- Lamp5_DESeq2_DEG$X
Astro_DESeq2_gene_list <- Astro_DESeq2_DEG$X
Peri_DESeq2_gene_list <- Peri_DESeq2_DEG$X
Endo_DESeq2_gene_list <- Endo_DESeq2_DEG$X
unique_DESeq2_genes <- unique(c(L2_3_IT_DESeq2_DEG$X,
                             L6_DESeq2_DEG$X,
                             Sst_DESeq2_DEG$X,
                             L5_DESeq2_DEG$X, 
                             L4_DESeq2_DEG$X,
                             Pvalb_DESeq2_DEG$X,
                             Sncg_DESeq2_DEG$X,
                             Non_neuronal_DESeq2_DEG$X,
                             Oligo_DESeq2_DEG$X,
                             Vip_DESeq2_DEG$X,
                             Lamp5_DESeq2_DEG$X,
                             Astro_DESeq2_DEG$X,
                             Peri_DESeq2_DEG$X,
                             Endo_DESeq2_DEG$X))

# List of genes differentially expressed per cluster for EdgeR
L2_3_IT_EdgeR_gene_list <- L2_3_IT_EdgeR_DEG$X
L6_EdgeR_gene_list <- L6_EdgeR_DEG$X
Sst_EdgeR_gene_list <- Sst_EdgeR_DEG$X
L5_EdgeR_gene_list <- L5_EdgeR_DEG$X
L4_EdgeR_gene_list <- L4_EdgeR_DEG$X
Pvalb_EdgeR_gene_list <- Pvalb_EdgeR_DEG$X
Sncg_EdgeR_gene_list <- Sncg_EdgeR_DEG$X
Non_neuronal_EdgeR_gene_list <- Non_neuronal_EdgeR_DEG$X
Oligo_EdgeR_gene_list <- Oligo_EdgeR_DEG$X
Vip_EdgeR_gene_list <- Vip_EdgeR_DEG$X
Lamp5_EdgeR_gene_list <- Lamp5_EdgeR_DEG$X
Astro_EdgeR_gene_list <- Astro_EdgeR_DEG$X
Peri_EdgeR_gene_list <- Peri_EdgeR_DEG$X
Endo_EdgeR_gene_list <- Endo_EdgeR_DEG$X
unique_EdgeR_genes <- unique(c(L2_3_IT_EdgeR_DEG$X,
                            L6_EdgeR_DEG$X,
                            Sst_EdgeR_DEG$X,
                            L5_EdgeR_DEG$X, 
                            L4_EdgeR_DEG$X,
                            Pvalb_EdgeR_DEG$X,
                            Sncg_EdgeR_DEG$X,
                            Non_neuronal_EdgeR_DEG$X,
                            Oligo_EdgeR_DEG$X,
                            Vip_EdgeR_DEG$X,
                            Lamp5_EdgeR_DEG$X,
                            Astro_EdgeR_DEG$X,
                            Peri_EdgeR_DEG$X,
                            Endo_EdgeR_DEG$X
                            ))

# Venn Diagram for Limma vs. DESeq2 vs. EdgeR per cluster
L2_3_IT_venn_list <- list(L2_3_IT_Limma_gene_list, L2_3_IT_EdgeR_gene_list, L2_3_IT_DESeq2_gene_list)
L2_3_IT_venn <-ggVennDiagram(L2_3_IT_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L2_3_IT", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L2_3_IT_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

L6_venn_list <- list(L6_Limma_gene_list, L6_EdgeR_gene_list, L6_DESeq2_gene_list)
L6_venn <- ggVennDiagram(L6_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L6", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L6_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Sst_venn_list <- list(Sst_Limma_gene_list, Sst_EdgeR_gene_list, Sst_DESeq2_gene_list)
Sst_venn <- ggVennDiagram(Sst_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Sst", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Sst_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

L5_venn_list <- list(L5_Limma_gene_list, L5_EdgeR_gene_list, L5_DESeq2_gene_list)
L5_venn <- ggVennDiagram(L5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L5", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L5_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

L4_venn_list <- list(L4_Limma_gene_list, L4_EdgeR_gene_list, L4_DESeq2_gene_list)
L4_venn <- ggVennDiagram(L4_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L4", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L4_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Pvalb_venn_list <- list(Pvalb_Limma_gene_list, Pvalb_EdgeR_gene_list, Pvalb_DESeq2_gene_list)
Pvalb_venn <- ggVennDiagram(Pvalb_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Pvalb", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Pvalb_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Sncg_venn_list <- list(Sncg_Limma_gene_list, Sncg_EdgeR_gene_list, Sncg_DESeq2_gene_list)
Sncg_venn <- ggVennDiagram(Sncg_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Sncg", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Sncg_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Non_neuronal_venn_list <- list(Non_neuronal_Limma_gene_list, Non_neuronal_EdgeR_gene_list, Non_neuronal_DESeq2_gene_list)
Non_neuronal_venn <- ggVennDiagram(Non_neuronal_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Non_neuronal", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Non_neuronal_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Oligo_venn_list <- list(Oligo_Limma_gene_list, Oligo_EdgeR_gene_list)
Oligo_venn <- ggVennDiagram(Oligo_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Oligo", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Oligo_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Vip_venn_list <- list(Vip_Limma_gene_list, Vip_EdgeR_gene_list, Vip_DESeq2_gene_list)
Vip_venn <- ggVennDiagram(Vip_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Vip", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Vip_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Lamp5_venn_list <- list(Lamp5_Limma_gene_list, Lamp5_EdgeR_gene_list, Lamp5_DESeq2_gene_list)
Lamp5_venn <- ggVennDiagram(Lamp5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Lamp5", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Lamp5_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Astro_venn_list <- list(Astro_Limma_gene_list, Astro_EdgeR_gene_list, Astro_DESeq2_gene_list)
Astro_venn <- ggVennDiagram(Astro_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Astro", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Astro_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Peri_venn_list <- list(Peri_Limma_gene_list, Peri_EdgeR_gene_list, Peri_DESeq2_gene_list)
Peri_venn <- ggVennDiagram(Peri_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Peri", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Peri_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

Endo_venn_list <- list(Endo_Limma_gene_list, Endo_EdgeR_gene_list, Endo_DESeq2_gene_list)
Endo_venn <- ggVennDiagram(Endo_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Endo", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Endo_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

unique_venn_list <- list(unique_Limma_genes, unique_EdgeR_genes, unique_DESeq2_genes)
unique_genes_venn <- ggVennDiagram(unique_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Unique Differentially Expressed Genes Identified for All Cell Types", subtitle = subtitle_info) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("unique_genes_M_MUT_and_WT_M_P60_CORT_venn.pdf", device = "pdf", path = venn_dir)

# Show genes identified by all methods for cell types
Reduce(intersect, list(L2_3_IT_Limma_gene_list, L2_3_IT_DESeq2_gene_list, L2_3_IT_EdgeR_gene_list))
Reduce(intersect, list(L6_Limma_gene_list, L6_DESeq2_gene_list, L6_EdgeR_gene_list))
Reduce(intersect, list(Sst_Limma_gene_list, Sst_DESeq2_gene_list, Sst_EdgeR_gene_list))
Reduce(intersect, list(L5_Limma_gene_list, L5_DESeq2_gene_list, L5_EdgeR_gene_list))
Reduce(intersect, list(L4_Limma_gene_list, L4_DESeq2_gene_list, L4_EdgeR_gene_list))
Reduce(intersect, list(Pvalb_Limma_gene_list, Pvalb_DESeq2_gene_list, Pvalb_EdgeR_gene_list))
Reduce(intersect, list(Sncg_Limma_gene_list, Sncg_DESeq2_gene_list, Sncg_EdgeR_gene_list))
Reduce(intersect, list(Non_neuronal_Limma_gene_list, Non_neuronal_DESeq2_gene_list, Non_neuronal_EdgeR_gene_list))
Reduce(intersect, list(Oligo_Limma_gene_list, Oligo_DESeq2_gene_list, Oligo_EdgeR_gene_list))
Reduce(intersect, list(Vip_Limma_gene_list, Vip_DESeq2_gene_list, Vip_EdgeR_gene_list))
Reduce(intersect, list(Lamp5_Limma_gene_list, Lamp5_DESeq2_gene_list, Lamp5_EdgeR_gene_list))
Reduce(intersect, list(Astro_Limma_gene_list, Astro_DESeq2_gene_list, Astro_EdgeR_gene_list))
Reduce(intersect, list(Peri_Limma_gene_list, Peri_DESeq2_gene_list, Peri_EdgeR_gene_list))
Reduce(intersect, list(Endo_Limma_gene_list, Endo_DESeq2_gene_list, Endo_EdgeR_gene_list))
