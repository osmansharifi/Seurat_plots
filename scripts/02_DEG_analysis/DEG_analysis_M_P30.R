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


################################################################################
# Venn Diagram for Differentially Expressed Genes Per Analysis
# By Viktoria Haghani

# Notes
# Sncg is excluded because no DEG were identified by any of the programs
# Oligo only has Limma and EdgeR since DESeq2 didn't identify any DEG
# Astro, Peri, and Endo only had DEG identified by one program each, so they are also excluded

# List of genes differentially expressed per cluster for Limma
L2_3_IT_Limma_gene_list <- L2_3_IT_Limma_stat_sig$X
L6_Limma_gene_list <- L6_Limma_stat_sig$X
Sst_Limma_gene_list <- Sst_Limma_stat_sig$X
L5_Limma_gene_list <- L5_Limma_stat_sig$X
L4_Limma_gene_list <- L4_Limma_stat_sig$X
Pvalb_Limma_gene_list <- Pvalb_Limma_stat_sig$X
Non_neuronal_Limma_gene_list <- Non_neuronal_Limma_stat_sig$X
Oligo_Limma_gene_list <- Oligo_Limma_stat_sig$X
Vip_Limma_gene_list <- Vip_Limma_stat_sig$X
Lamp5_Limma_gene_list <- Lamp5_Limma_stat_sig$X
Astro_Limma_gene_list <- Astro_Limma_stat_sig$X
Peri_Limma_gene_list <- Peri_Limma_stat_sig$X
Endo_Limma_gene_list <- Endo_Limma_stat_sig$X
unique_Limma_genes <- unique(c(L2_3_IT_Limma_stat_sig$X,
                     L6_Limma_stat_sig$X,
                     Sst_Limma_stat_sig$X,
                     L5_Limma_stat_sig$X, 
                     L4_Limma_stat_sig$X,
                     Pvalb_Limma_stat_sig$X,
                     Non_neuronal_Limma_stat_sig$X,
                     Oligo_Limma_stat_sig$X,
                     Vip_Limma_stat_sig$X,
                     Lamp5_Limma_stat_sig$X,
                     Astro_Limma_stat_sig$X,
                     Peri_Limma_stat_sig$X,
                     Endo_Limma_stat_sig$X))
all_Limma_genes_not_unique <- c(L2_3_IT_Limma_stat_sig$X,
                            L6_Limma_stat_sig$X,
                            Sst_Limma_stat_sig$X,
                            L5_Limma_stat_sig$X, 
                            L4_Limma_stat_sig$X,
                            Pvalb_Limma_stat_sig$X,
                            Non_neuronal_Limma_stat_sig$X,
                            Oligo_Limma_stat_sig$X,
                            Vip_Limma_stat_sig$X,
                            Lamp5_Limma_stat_sig$X,
                            Astro_Limma_stat_sig$X,
                            Peri_Limma_stat_sig$X,
                            Endo_Limma_stat_sig$X)

# List of genes differentially expressed per cluster for DESeq2
L2_3_IT_DESeq2_gene_list <- L2_3_IT_DESeq2_DEG_stat_sig$X
L6_DESeq2_gene_list <- L6_DESeq2_DEG_stat_sig$X
Sst_DESeq2_gene_list <- Sst_DESeq2_DEG_stat_sig$X
L5_DESeq2_gene_list <- L5_DESeq2_DEG_stat_sig$X
L4_DESeq2_gene_list <- L4_DESeq2_DEG_stat_sig$X
Pvalb_DESeq2_gene_list <- Pvalb_DESeq2_DEG_stat_sig$X
Non_neuronal_DESeq2_gene_list <- Non_neuronal_DESeq2_DEG_stat_sig$X
Oligo_DESeq2_gene_list <- Oligo_DESeq2_DEG_stat_sig$X
Vip_DESeq2_gene_list <- Vip_DESeq2_DEG_stat_sig$X
Lamp5_DESeq2_gene_list <- Lamp5_DESeq2_DEG_stat_sig$X
Astro_DESeq2_gene_list <- Astro_DESeq2_DEG_stat_sig$X
Peri_DESeq2_gene_list <- Peri_DESeq2_DEG_stat_sig$X
Endo_DESeq2_gene_list <- Endo_DESeq2_DEG_stat_sig$X
unique_DESeq2_genes <- unique(c(L2_3_IT_DESeq2_DEG_stat_sig$X,
                             L6_DESeq2_DEG_stat_sig$X,
                             Sst_DESeq2_DEG_stat_sig$X,
                             L5_DESeq2_DEG_stat_sig$X, 
                             L4_DESeq2_DEG_stat_sig$X,
                             Pvalb_DESeq2_DEG_stat_sig$X,
                             Non_neuronal_DESeq2_DEG_stat_sig$X,
                             Oligo_DESeq2_DEG_stat_sig$X,
                             Vip_DESeq2_DEG_stat_sig$X,
                             Lamp5_DESeq2_DEG_stat_sig$X,
                             Astro_DESeq2_DEG_stat_sig$X,
                             Peri_DESeq2_DEG_stat_sig$X,
                             Endo_DESeq2_DEG_stat_sig$X))

# List of genes differentially expressed per cluster for EdgeR
L2_3_IT_EdgeR_gene_list <- L2_3_IT_EdgeR_stat_sig$X
L6_EdgeR_gene_list <- L6_EdgeR_stat_sig$X
Sst_EdgeR_gene_list <- Sst_EdgeR_stat_sig$X
L5_EdgeR_gene_list <- L5_EdgeR_stat_sig$X
L4_EdgeR_gene_list <- L4_EdgeR_stat_sig$X
Pvalb_EdgeR_gene_list <- Pvalb_EdgeR_stat_sig$X
Non_neuronal_EdgeR_gene_list <- Non_neuronal_EdgeR_stat_sig$X
Oligo_EdgeR_gene_list <- Oligo_EdgeR_stat_sig$X
Vip_EdgeR_gene_list <- Vip_EdgeR_stat_sig$X
Lamp5_EdgeR_gene_list <- Lamp5_EdgeR_stat_sig$X
Astro_EdgeR_gene_list <- Astro_EdgeR_stat_sig$X
Peri_EdgeR_gene_list <- Peri_EdgeR_stat_sig$X
Endo_EdgeR_gene_list <- Endo_EdgeR_stat_sig$X
unique_EdgeR_genes <- unique(c(L2_3_IT_EdgeR_stat_sig$X,
                            L6_EdgeR_stat_sig$X,
                            Sst_EdgeR_stat_sig$X,
                            L5_EdgeR_stat_sig$X, 
                            L4_EdgeR_stat_sig$X,
                            Pvalb_EdgeR_stat_sig$X,
                            Non_neuronal_EdgeR_stat_sig$X,
                            Oligo_EdgeR_stat_sig$X,
                            Vip_EdgeR_stat_sig$X,
                            Lamp5_EdgeR_stat_sig$X,
                            Astro_EdgeR_stat_sig$X,
                            Peri_EdgeR_stat_sig$X,
                            Endo_EdgeR_stat_sig$X
                            ))

# Venn Diagram for Limma vs. DESeq2 vs. EdgeR per cluster
L2_3_IT_venn_list <- list(L2_3_IT_Limma_gene_list, L2_3_IT_EdgeR_gene_list, L2_3_IT_DESeq2_gene_list)
#L2_3_IT_venn <- 
ggVennDiagram(L2_3_IT_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L2_3_IT") +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave("L2_3_IT_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

L6_venn_list <- list(L6_Limma_gene_list, L6_EdgeR_gene_list, L6_DESeq2_gene_list)
L6_venn <- ggVennDiagram(L6_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L6") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L6_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Sst_venn_list <- list(Sst_Limma_gene_list, Sst_EdgeR_gene_list, Sst_DESeq2_gene_list)
Sst_venn <- ggVennDiagram(Sst_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Sst") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Sst_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

L5_venn_list <- list(L5_Limma_gene_list, L5_EdgeR_gene_list, L5_DESeq2_gene_list)
L5_venn <- ggVennDiagram(L5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L5_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

L4_venn_list <- list(L4_Limma_gene_list, L4_EdgeR_gene_list, L4_DESeq2_gene_list)
L4_venn <- ggVennDiagram(L4_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for L4") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("L4_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Pvalb_venn_list <- list(Pvalb_Limma_gene_list, Pvalb_EdgeR_gene_list, Pvalb_DESeq2_gene_list)
Pvalb_venn <- ggVennDiagram(Pvalb_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Pvalb") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Pvalb_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Non_neuronal_venn_list <- list(Non_neuronal_Limma_gene_list, Non_neuronal_EdgeR_gene_list, Non_neuronal_DESeq2_gene_list)
Non_neuronal_venn <- ggVennDiagram(Non_neuronal_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Non_neuronal") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Non_neuronal_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Oligo_venn_list <- list(Oligo_Limma_gene_list, Oligo_EdgeR_gene_list)
Oligo_venn <- ggVennDiagram(Oligo_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Oligo") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Oligo_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Vip_venn_list <- list(Vip_Limma_gene_list, Vip_EdgeR_gene_list, Vip_DESeq2_gene_list)
Vip_venn <- ggVennDiagram(Vip_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Vip") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Vip_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

Lamp5_venn_list <- list(Lamp5_Limma_gene_list, Lamp5_EdgeR_gene_list, Lamp5_DESeq2_gene_list)
Lamp5_venn <- ggVennDiagram(Lamp5_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Differentially Expressed Genes Identified for Lamp5") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Lamp5_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

unique_venn_list <- list(unique_Limma_genes, unique_EdgeR_genes, unique_DESeq2_genes)
unique_genes_venn <- ggVennDiagram(unique_venn_list, color = "black", lwd = 0.8, lty = 1, category.names = c("Limma", "EdgeR", "DESeq2")) +
  ggplot2::scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Unique Differentially Expressed Genes Identified for All Cell Types") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("unique_genes_venn.pdf", device = "pdf", path = "~/GitHub/snRNA-seq-pipeline/DEG_data/venn_diagrams")

# Show genes identified by all methods for cell types
Reduce(intersect, list(L2_3_IT_Limma_gene_list, L2_3_IT_DESeq2_gene_list, L2_3_IT_EdgeR_gene_list))
Reduce(intersect, list(L6_Limma_gene_list, L6_DESeq2_gene_list, L6_EdgeR_gene_list))
Reduce(intersect, list(Sst_Limma_gene_list, Sst_DESeq2_gene_list, Sst_EdgeR_gene_list))
Reduce(intersect, list(L5_Limma_gene_list, L5_DESeq2_gene_list, L5_EdgeR_gene_list))
Reduce(intersect, list(L4_Limma_gene_list, L4_DESeq2_gene_list, L4_EdgeR_gene_list))
Reduce(intersect, list(Pvalb_Limma_gene_list, Pvalb_DESeq2_gene_list, Pvalb_EdgeR_gene_list))
Reduce(intersect, list(Vip_Limma_gene_list, Vip_DESeq2_gene_list, Vip_EdgeR_gene_list))
Reduce(intersect, list(Lamp5_Limma_gene_list, Lamp5_DESeq2_gene_list, Lamp5_EdgeR_gene_list))
