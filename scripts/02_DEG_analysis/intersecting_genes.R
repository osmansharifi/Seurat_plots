library(glue)

# Paths
EdgeR_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/EdgeR/M_MUT_and_WT_M_E18_WB/"
DESeq2_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/DESeq2/M_MUT_and_WT_M_E18_WB/"
Limma_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_E18_WB/"
DEG_data_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/all_methods/M_MUT_and_WT_M_E18_WB/"

# Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 

# Other variables
metadata_info <- "M_MUT_and_WT_M_E18_WB"
subtitle_info <- "Mice, Male, E18, Whole Brain"

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

# Find genes identified by all methods for cell types
L2_3_IT_all_methods <- Reduce(intersect, list(L2_3_IT_Limma_gene_list, L2_3_IT_DESeq2_gene_list, L2_3_IT_EdgeR_gene_list))
L2_3_IT_all_methods_data <- L2_3_IT_EdgeR_DEG[L2_3_IT_EdgeR_DEG$X %in% L2_3_IT_all_methods, ]
write.csv(L2_3_IT_all_methods_data, file = glue(DEG_data_dir, "L2_3_IT", "_", metadata_info, "_all_methods_DEG.csv"))

L6_all_methods <- Reduce(intersect, list(L6_Limma_gene_list, L6_DESeq2_gene_list, L6_EdgeR_gene_list))
L6_all_methods_data <- L6_EdgeR_DEG[L6_EdgeR_DEG$X %in% L6_all_methods, ]
write.csv(L6_all_methods_data, file = glue(DEG_data_dir, "L6", "_", metadata_info, "_all_methods_DEG.csv"))

Sst_all_methods <- Reduce(intersect, list(Sst_Limma_gene_list, Sst_DESeq2_gene_list, Sst_EdgeR_gene_list))
Sst_all_methods_data <- Sst_EdgeR_DEG[Sst_EdgeR_DEG$X %in% Sst_all_methods, ]
write.csv(Sst_all_methods_data, file = glue(DEG_data_dir, "Sst", "_", metadata_info, "_all_methods_DEG.csv"))

L5_all_methods <- Reduce(intersect, list(L5_Limma_gene_list, L5_DESeq2_gene_list, L5_EdgeR_gene_list))
L5_all_methods_data <- L5_EdgeR_DEG[L5_EdgeR_DEG$X %in% L5_all_methods, ]
write.csv(L5_all_methods_data, file = glue(DEG_data_dir, "L5", "_", metadata_info, "_all_methods_DEG.csv"))

L4_all_methods <- Reduce(intersect, list(L4_Limma_gene_list, L4_DESeq2_gene_list, L4_EdgeR_gene_list))
L4_all_methods_data <- L4_EdgeR_DEG[L4_EdgeR_DEG$X %in% L4_all_methods, ]
write.csv(L4_all_methods_data, file = glue(DEG_data_dir, "L4", "_", metadata_info, "_all_methods_DEG.csv"))

Pvalb_all_methods <- Reduce(intersect, list(Pvalb_Limma_gene_list, Pvalb_DESeq2_gene_list, Pvalb_EdgeR_gene_list))
Pvalb_all_methods_data <- Pvalb_EdgeR_DEG[Pvalb_EdgeR_DEG$X %in% Pvalb_all_methods, ]
write.csv(Pvalb_all_methods_data, file = glue(DEG_data_dir, "Pvalb", "_", metadata_info, "_all_methods_DEG.csv"))

Sncg_all_methods <- Reduce(intersect, list(Sncg_Limma_gene_list, Sncg_DESeq2_gene_list, Sncg_EdgeR_gene_list))
Sncg_all_methods_data <- Sncg_EdgeR_DEG[Sncg_EdgeR_DEG$X %in% Sncg_all_methods, ]
write.csv(Sncg_all_methods_data, file = glue(DEG_data_dir, "Sncg", "_", metadata_info, "_all_methods_DEG.csv"))

Non_neuronal_all_methods <- Reduce(intersect, list(Non_neuronal_Limma_gene_list, Non_neuronal_DESeq2_gene_list, Non_neuronal_EdgeR_gene_list))
Non_neuronal_all_methods_data <- Non_neuronal_EdgeR_DEG[Non_neuronal_EdgeR_DEG$X %in% Non_neuronal_all_methods, ]
write.csv(Non_neuronal_all_methods_data, file = glue(DEG_data_dir, "Non_neuronal", "_", metadata_info, "_all_methods_DEG.csv"))

Oligo_all_methods <- Reduce(intersect, list(Oligo_Limma_gene_list, Oligo_DESeq2_gene_list, Oligo_EdgeR_gene_list))
Oligo_all_methods_data <- Oligo_EdgeR_DEG[Oligo_EdgeR_DEG$X %in% Oligo_all_methods, ]
write.csv(Oligo_all_methods_data, file = glue(DEG_data_dir, "Oligo", "_", metadata_info, "_all_methods_DEG.csv"))

Vip_all_methods <- Reduce(intersect, list(Vip_Limma_gene_list, Vip_DESeq2_gene_list, Vip_EdgeR_gene_list))
Vip_all_methods_data <- Vip_EdgeR_DEG[Vip_EdgeR_DEG$X %in% Vip_all_methods, ]
write.csv(Vip_all_methods_data, file = glue(DEG_data_dir, "Vip", "_", metadata_info, "_all_methods_DEG.csv"))

Lamp5_all_methods <- Reduce(intersect, list(Lamp5_Limma_gene_list, Lamp5_DESeq2_gene_list, Lamp5_EdgeR_gene_list))
Lamp5_all_methods_data <- Lamp5_EdgeR_DEG[Lamp5_EdgeR_DEG$X %in% Lamp5_all_methods, ]
write.csv(Lamp5_all_methods_data, file = glue(DEG_data_dir, "Lamp5", "_", metadata_info, "_all_methods_DEG.csv"))

Astro_all_methods <- Reduce(intersect, list(Astro_Limma_gene_list, Astro_DESeq2_gene_list, Astro_EdgeR_gene_list))
Astro_all_methods_data <- Astro_EdgeR_DEG[Astro_EdgeR_DEG$X %in% Astro_all_methods, ]
write.csv(Astro_all_methods_data, file = glue(DEG_data_dir, "Astro", "_", metadata_info, "_all_methods_DEG.csv"))

Peri_all_methods <- Reduce(intersect, list(Peri_Limma_gene_list, Peri_DESeq2_gene_list, Peri_EdgeR_gene_list))
Peri_all_methods_data <- Peri_EdgeR_DEG[Peri_EdgeR_DEG$X %in% Peri_all_methods, ]
write.csv(Peri_all_methods_data, file = glue(DEG_data_dir, "Peri", "_", metadata_info, "_all_methods_DEG.csv"))

Endo_all_methods <- Reduce(intersect, list(Endo_Limma_gene_list, Endo_DESeq2_gene_list, Endo_EdgeR_gene_list))
Endo_all_methods_data <- Endo_EdgeR_DEG[Endo_EdgeR_DEG$X %in% Endo_all_methods, ]
write.csv(Endo_all_methods_data, file = glue(DEG_data_dir, "Endo", "_", metadata_info, "_all_methods_DEG.csv"))