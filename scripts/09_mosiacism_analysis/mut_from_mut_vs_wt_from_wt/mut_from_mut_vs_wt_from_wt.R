########################################################
## Create broad categories and run mosiacism analysis ##
########################################################
library(magrittr)
library(VennDiagram)
library(grDevices)
library(dplyr)
library(Seurat)
library(glue)
library(scCustomize)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(GeneOverlap)
library(tidyverse)

##################
## Load samples ##
##################
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
load(glue('{base_path}/all.female.cortex.parsed.RData'))
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/mut_from_mut_vs_wt_from_wt'

##################################################################
## select the mut cells from het females and wt from wt females ##
##################################################################
# Create a new column called broad_class based on celltype.call
all.female.cortex$broad_class <- ifelse(
  all.female.cortex$celltype.call %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Vip"), 
  "GABAergic", 
  ifelse(
    all.female.cortex$celltype.call %in% c("L2_3_IT", "L4", "L5", "L6"), 
    "Glutamatergic", 
    ifelse(
      all.female.cortex$celltype.call %in% c("Astro", "Non-neuronal", "Oligo"), 
      "Non-neuronal", 
      "Other"
    )
  )
)
MUT_from_MUT = Cells(all.female.cortex)[which(all.female.cortex$Condition == "MUTANT" & all.female.cortex$Mecp2_allele == "Mecp2_MUT" & all.female.cortex$broad_class == "GABAergic")]
WT_from_WT = Cells(all.female.cortex)[which(all.female.cortex$Condition == "WT" & all.female.cortex$Mecp2_allele == "Mecp2_WT" & all.female.cortex$broad_class == "GABAergic")]
slct_MUT_from_MUT = sample(MUT_from_MUT, size = 135)
slct_WT_from_WT = sample(WT_from_WT, size = 135)
subset_cell_nonautonomous = subset(all.female.cortex, cells = c(slct_MUT_from_MUT, slct_WT_from_WT))

DimPlot_scCustom(seurat_object = subset_cell_nonautonomous, group.by = "Mecp2_allele", split.by = 'broad_class', pt.size = 0.8)
ggplot2::ggsave(glue("{base_path}/celltype_Mecp2_allele_umap.pdf"),
                device = NULL,
                height = 8.5,
                width = 12,
                plot = last_plot())
#########################################################################################################
## Perform DEG analysis between the MUT cells from the HET mouse and WT cells from the WT mouse cortex ##
#########################################################################################################
age_groups <- unique(subset_cell_nonautonomous@meta.data$Age)
deg_results <- list()

for (age_group in age_groups) {
  cat("Performing DEG analysis for", age_group, "\n")
  
  # Subset cells based on age
  age_subset <- subset(subset_cell_nonautonomous, subset = Age == age_group)
  
  # Get expression info
  expr <- as.matrix(GetAssayData(age_subset))
  
  # Filter out genes that are 0 for every cell
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad, ]
  
  logcpm <- cpm(expr, prior.count = 2, log = TRUE)
  mm <- model.matrix(~0 + Condition, data = age_subset@meta.data)
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Extract DEG results
  contrasts <- makeContrasts(c(ConditionMUTANT) - c(ConditionWT), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contrasts)
  tmp <- eBayes(tmp)
  top_table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
  
  # Store DEG results for this age group
  deg_results[[age_group]] <- top_table
}

deg_results$P30$Timepoint <- 'P30'
deg_results$P60$Timepoint <- 'P60'
deg_results$P150$Timepoint <- 'P150'
deg_results$P30$DEG_Test <- 'MUT cells from HET vs WT cells from WT'
deg_results$P60$DEG_Test <- 'MUT cells from HET vs WT cells from WT'
deg_results$P150$DEG_Test <- 'MUT cells from HET vs WT cells from WT'
MUTvsWT <- rbind(deg_results$P30, deg_results$P60, deg_results$P150)
MUTvsWT$SYMBOL <- rownames(MUTvsWT)
write.csv(MUTvsWT, file = glue('{base_path}/mutvswt_DEG_GABA.csv'))
# Create an empty data frame to store results
result_df <- data.frame()

for (age_group in age_groups) {
  cat("DEG analysis results for", age_group, "\n")
  
  deg_table <- deg_results[[age_group]]
  
  # Add a new column "Time_point" with the current age_group
  deg_table$Time_point <- age_group
  
  # Filter rows where adj.P.Val < 0.05
  deg_table_filtered <- deg_table[deg_table$adj.P.Val < 0.05, ]
  
  num_degs <- nrow(deg_table_filtered)
  print(num_degs)
  
  # Append the filtered deg_table to the result_df
  # Now result_df contains only rows with adj.P.Val < 0.05
  result_df <- rbind(result_df, deg_table_filtered)
}
####################
## Write csv file ##
####################
# Assuming result_df is not empty
if (nrow(result_df) > 0) {
  # Specify the file path where you want to save the CSV file
  csv_file_path <- glue("{base_path}/sig_DEGs_MUTvsWT_gaba.csv")
  
  # Write result_df to CSV
  write.csv(result_df, file = csv_file_path, row.names = FALSE)
  
  cat("result_df has been successfully written to", csv_file_path, "\n")
} else {
  cat("result_df is empty, nothing to write to CSV.\n")
}

top.table <- deg_results$P60
top.table$Gene <- rownames(top.table)
top.table$diffexpressed <- 'NO'
top.table$diffexpressed[top.table$logFC > 0 & top.table$adj.P.Val < 0.05] <- 'UP'
top.table$diffexpressed[top.table$logFC < 0 & top.table$adj.P.Val < 0.05] <- 'DOWN'
top.table$diffexpressed[top.table$adj.P.Val > 0.05] <- 'Not Sig'
top.table$delabel <- NA
thresh = head(arrange(top.table, adj.P.Val), 10)$adj.P.Val[10]
top.table$delabel[top.table$adj.P.Val <=thresh] <-(top.table$Gene[top.table$adj.P.Val<=thresh])

# Volcano Plot
ggplot(data = top.table, aes(x = logFC, y = -log2(adj.P.Val), col = diffexpressed, label = delabel)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf, box.padding = 0.8) +
  scale_color_manual(values = c('blue', 'black', 'red')) +
  #scale_color_manual(values = c('black')) +
  scale_y_continuous(limits = c(0, 18)) +  # Set the limits here
  theme(
    text = element_text(size = 16),
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14, face = "bold"),
    title = element_text(size = 14, face = "bold")
  ) +
  labs(title = 'P30 GABAergic MUT cells from MUT females vs WT cells from WT female')
ggplot2::ggsave(glue("{base_path}/Vol_gaba_MUTvsWT_P30females.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

###########################################
## Venn diagram of the overlapping genes ##
###########################################

# Extract significant DEGs
sig_genes_P30 <- rownames(deg_results$P30)[deg_results$P30$adj.P.Val < 0.05]
sig_genes_P60 <- rownames(deg_results$P60)[deg_results$P60$adj.P.Val < 0.05]
sig_genes_P150 <- rownames(deg_results$P150)[deg_results$P150$adj.P.Val < 0.05]
intersection_all2 <- intersect(sig_genes_P150,sig_genes_P30)
intersection_all3 <- intersect(intersection_all2,sig_genes_P60)
intersection_all4 <- intersect(sig_genes_P30,sig_genes_P60)
# Create a Venn diagram
pdf(glue("{base_path}/venn_GABAergic.pdf"))
temp <- venn.diagram(
  x = list(
    P30 = sig_genes_P30,
    P60 = sig_genes_P60,
    P150 = sig_genes_P150
  ),
  category.names = c("P30", "P60", "P150"),
  main = 'GABAergic DEGs from MUT cells from MUT females and WT cells from WT females ',
  #filename = glue("{base_path}/broad_group_analysis/venn_glutamatergic.pdf"),
  filename = NULL,
  col = c('#E6B8BFFF', '#CC7A88FF', '#990F26FF'), 
  fill = c('#E6B8BFFF', '#CC7A88FF', '#990F26FF'),
  cat.cex = 1.2,
  cat.fontface = "bold",
  euler.d = TRUE,
  disable.logging = TRUE,
  hyper.test = TRUE
)
grid.draw(temp)
dev.off()

# Make all vectors the same length (pad with NA)
max_length <- max(length(sig_genes_P30), length(sig_genes_P60), length(sig_genes_P150))
sig_genes_P30 <- c(sig_genes_P30, rep(NA, max_length - length(sig_genes_P30)))
sig_genes_P60 <- c(sig_genes_P60, rep(NA, max_length - length(sig_genes_P60)))
sig_genes_P150 <- c(sig_genes_P150, rep(NA, max_length - length(sig_genes_P150)))

# Combine into a dataframe
combined_df <- data.frame(sig_genes_P30, sig_genes_P60, sig_genes_P150)

go.obj <- newGeneOverlap(combined_df$sig_genes_P30,
                         combined_df$sig_genes_P60,
                         combined_df$sig_genes_P150,
                         genome.size = 22000)
go.obj <- testGeneOverlap(go.obj)
getPval(go.obj)
getOddsRatio(go.obj)
getJaccard(go.obj)
getContbl(go.obj)
print(go.obj)
write.csv(print(go.obj), file = glue('{base_path}/geneoverlap_gaba_3timepoints.txt'))
DEGs = top.table
# Top DEGs

DEGs <- deg_results$P60 %>%
  tibble::rownames_to_column() %>%
  tibble::as_tibble() %>%
  dplyr::rename(SYMBOL = rowname) %>%
  dplyr::mutate(FC = dplyr::case_when(logFC >0 ~ 2^logFC,
                                      logFC <0 ~ -1/(2^logFC))) %>%
  dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %T>%
  openxlsx::write.xlsx(file=glue::glue("{base_path}/DEGs.xlsx")) %>%
  dplyr::filter(P.Value < 0.05) %T>%
  openxlsx::write.xlsx(file=glue::glue("{base_path}/sig_DEGs.xlsx"))

enrichR:::.onAttach()
source(glue('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/GO_ploting_functions.R'))
tryCatch({
  DEGs %>% 
    dplyr::select(SYMBOL) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Molecular_Function_2018",
                       "GO_Cellular_Component_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
    openxlsx::write.xlsx(file=glue::glue("{base_path}/GABA_P60_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Mm.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{base_path}/GABA_P60_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("{base_path}/GABA_P60_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
  })
print(glue::glue("The pipeline has finished for samples"))
