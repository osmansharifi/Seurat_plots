library(foreach)
library(glue)
library(ggplot2)

# By Osman Sharifi & Viktoria Haghani

# Note that GO_analysis_enrichment_scores.R should be run before running this script
# This script reads in the GO Terms previously identified

################################################################################
## Variables

#Paths
go_data_path <- "~/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/"
figure_path <- "~/GitHub/snRNA-seq-pipeline/figures/go_analysis/enrichment_scores/"

# Lists
time_points <- list("E18_WB", "P30_CORT", "P60_CORT", "P120_CORT")
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo")
go_ontologies <- list("BP", "CC", "MF")

# Variable
metadata <- "M_MUT_and_WT_M"
plot_title <- "GO Analysis"
plot_subtitle <- "All Cell Types, Time Points, and Ontologies for Male Mice"

################################################################################

# Read in GO Data
for (time_point in time_points){
  for (cell_type in cell_types){
    for (ont in go_ontologies){
      # Read in all GO Data
      assign(glue(cell_type, "_", time_point, "_", ont), read.csv(file = glue(go_data_path, metadata, "_", time_point, "/", cell_type, "_", metadata, "_", time_point, "_", ont, "_gentable.csv")))
    }
  }
}

go_data <- list(Astro_E18_WB_BP,
            Astro_E18_WB_CC,
            Astro_E18_WB_MF,
            Endo_E18_WB_BP,
            Endo_E18_WB_CC,
            Endo_E18_WB_MF,
            L2_3_IT_E18_WB_BP,
            L2_3_IT_E18_WB_CC,
            L2_3_IT_E18_WB_MF,
            L4_E18_WB_BP,
            L4_E18_WB_CC,
            L4_E18_WB_MF,
            L5_E18_WB_BP,
            L5_E18_WB_CC,
            L5_E18_WB_MF,
            L6_E18_WB_BP,
            L6_E18_WB_CC,
            L6_E18_WB_MF,
            Lamp5_E18_WB_BP,
            Lamp5_E18_WB_CC,
            Lamp5_E18_WB_MF,
            Non_neuronal_E18_WB_BP,
            Non_neuronal_E18_WB_CC,
            Non_neuronal_E18_WB_MF,
            Oligo_E18_WB_BP,
            Oligo_E18_WB_CC,
            Oligo_E18_WB_MF,
            Peri_E18_WB_BP,
            Peri_E18_WB_CC,
            Peri_E18_WB_MF,
            Pvalb_E18_WB_BP,
            Pvalb_E18_WB_CC,
            Pvalb_E18_WB_MF,
            Sncg_E18_WB_BP,
            Sncg_E18_WB_CC,
            Sncg_E18_WB_MF,
            Sst_E18_WB_BP,
            Sst_E18_WB_CC,
            Sst_E18_WB_MF,
            Vip_E18_WB_BP,
            Vip_E18_WB_CC,
            Vip_E18_WB_MF,
            Astro_P30_CORT_BP,
            Astro_P30_CORT_CC,
            Astro_P30_CORT_MF,
            Endo_P30_CORT_BP,
            Endo_P30_CORT_CC,
            Endo_P30_CORT_MF,
            L2_3_IT_P30_CORT_BP,
            L2_3_IT_P30_CORT_CC,
            L2_3_IT_P30_CORT_MF,
            L4_P30_CORT_BP,
            L4_P30_CORT_CC,
            L4_P30_CORT_MF,
            L5_P30_CORT_BP,
            L5_P30_CORT_CC,
            L5_P30_CORT_MF,
            L6_P30_CORT_BP,
            L6_P30_CORT_CC,
            L6_P30_CORT_MF,
            Lamp5_P30_CORT_BP,
            Lamp5_P30_CORT_CC,
            Lamp5_P30_CORT_MF,
            Non_neuronal_P30_CORT_BP,
            Non_neuronal_P30_CORT_CC,
            Non_neuronal_P30_CORT_MF,
            Oligo_P30_CORT_BP,
            Oligo_P30_CORT_CC,
            Oligo_P30_CORT_MF,
            Peri_P30_CORT_BP,
            Peri_P30_CORT_CC,
            Peri_P30_CORT_MF,
            Pvalb_P30_CORT_BP,
            Pvalb_P30_CORT_CC,
            Pvalb_P30_CORT_MF,
            Sncg_P30_CORT_BP,
            Sncg_P30_CORT_CC,
            Sncg_P30_CORT_MF,
            Sst_P30_CORT_BP,
            Sst_P30_CORT_CC,
            Sst_P30_CORT_MF,
            Vip_P30_CORT_BP,
            Vip_P30_CORT_CC,
            Vip_P30_CORT_MF,
            Astro_P60_CORT_BP,
            Astro_P60_CORT_CC,
            Astro_P60_CORT_MF,
            Endo_P60_CORT_BP,
            Endo_P60_CORT_CC,
            Endo_P60_CORT_MF,
            L2_3_IT_P60_CORT_BP,
            L2_3_IT_P60_CORT_CC,
            L2_3_IT_P60_CORT_MF,
            L4_P60_CORT_BP,
            L4_P60_CORT_CC,
            L4_P60_CORT_MF,
            L5_P60_CORT_BP,
            L5_P60_CORT_CC,
            L5_P60_CORT_MF,
            L6_P60_CORT_BP,
            L6_P60_CORT_CC,
            L6_P60_CORT_MF,
            Lamp5_P60_CORT_BP,
            Lamp5_P60_CORT_CC,
            Lamp5_P60_CORT_MF,
            Non_neuronal_P60_CORT_BP,
            Non_neuronal_P60_CORT_CC,
            Non_neuronal_P60_CORT_MF,
            Oligo_P60_CORT_BP,
            Oligo_P60_CORT_CC,
            Oligo_P60_CORT_MF,
            Peri_P60_CORT_BP,
            Peri_P60_CORT_CC,
            Peri_P60_CORT_MF,
            Pvalb_P60_CORT_BP,
            Pvalb_P60_CORT_CC,
            Pvalb_P60_CORT_MF,
            Sncg_P60_CORT_BP,
            Sncg_P60_CORT_CC,
            Sncg_P60_CORT_MF,
            Sst_P60_CORT_BP,
            Sst_P60_CORT_CC,
            Sst_P60_CORT_MF,
            Vip_P60_CORT_BP,
            Vip_P60_CORT_CC,
            Vip_P60_CORT_MF,
            Astro_P120_CORT_BP,
            Astro_P120_CORT_CC,
            Astro_P120_CORT_MF,
            Endo_P120_CORT_BP,
            Endo_P120_CORT_CC,
            Endo_P120_CORT_MF,
            L2_3_IT_P120_CORT_BP,
            L2_3_IT_P120_CORT_CC,
            L2_3_IT_P120_CORT_MF,
            L4_P120_CORT_BP,
            L4_P120_CORT_CC,
            L4_P120_CORT_MF,
            L5_P120_CORT_BP,
            L5_P120_CORT_CC,
            L5_P120_CORT_MF,
            L6_P120_CORT_BP,
            L6_P120_CORT_CC,
            L6_P120_CORT_MF,
            Lamp5_P120_CORT_BP,
            Lamp5_P120_CORT_CC,
            Lamp5_P120_CORT_MF,
            Non_neuronal_P120_CORT_BP,
            Non_neuronal_P120_CORT_CC,
            Non_neuronal_P120_CORT_MF,
            Oligo_P120_CORT_BP,
            Oligo_P120_CORT_CC,
            Oligo_P120_CORT_MF,
            Peri_P120_CORT_BP,
            Peri_P120_CORT_CC,
            Peri_P120_CORT_MF,
            Pvalb_P120_CORT_BP,
            Pvalb_P120_CORT_CC,
            Pvalb_P120_CORT_MF,
            Sncg_P120_CORT_BP,
            Sncg_P120_CORT_CC,
            Sncg_P120_CORT_MF,
            Sst_P120_CORT_BP,
            Sst_P120_CORT_CC,
            Sst_P120_CORT_MF,
            Vip_P120_CORT_BP,
            Vip_P120_CORT_CC,
            Vip_P120_CORT_MF)

# create fake data
set.seed(1024) # keep reproducibility
go <- paste0("GO", sample(1000:2000, 5))
data <- data.frame("GOs" = rep(go, 2), 
                   "Condition" = rep(c("A", "B"), each = 5),
                   "GeneRatio" = 1 / sample(10, 10), 
                   "p.adjust" = 0.05 / sample(10, 10))

# plot: dot plot
ggplot(data = data, aes(x = Condition, y = GOs, 
                        color = `p.adjust`, size = GeneRatio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("Cell Types for Each Time Point") + 
  labs(title = plot_title,
       subtitle = plot_subtitle)