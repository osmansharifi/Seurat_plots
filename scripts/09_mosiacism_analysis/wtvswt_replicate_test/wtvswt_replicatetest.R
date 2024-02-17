###########################################
## Run replicate test between WT females ##
###########################################
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
load('/Users/osman/Desktop/LaSalle_lab/Seurat_objects/postnatal_cortex_20221109.RData')
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/wtvswt_replicate_test'
DefaultAssay(postnatal_cortex) <- "integrated"
##################################################################
## select the wt cells from wt females and wt from wt females ##
##################################################################
# filter samples 
postnatal_cortex <- subset(postnatal_cortex, subset = Age == 'P150' & Condition == 'WT')
postnatal_cortex$broad_class <- ifelse(
  postnatal_cortex$celltype.call %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Vip"), 
  "GABAergic", 
  ifelse(
    postnatal_cortex$celltype.call %in% c("L2_3_IT", "L4", "L5", "L6"), 
    "Glutamatergic", 
    ifelse(
      postnatal_cortex$celltype.call %in% c("Astro", "Non-neuronal", "Oligo"), 
      "Non-neuronal", 
      "Other"
    )
  )
)

# Create a new column in the metadata called "replicates"
postnatal_cortex$replicates <- NA

# Assign values based on conditions
postnatal_cortex$replicates[postnatal_cortex$orig.ident %in% c("WT_F_P150_CORT1", "WT_F_P150_CORT3")] <- "rep1rep2"
postnatal_cortex$replicates[postnatal_cortex$orig.ident %in% c("WT_F_P150_CORT2", "WT_F_P150_CORT4")] <- "rep3rep4"
# Select cells
WT1_WT3 = Cells(postnatal_cortex)[which(postnatal_cortex$replicates == "rep1rep2")]
WT2_WT4 = Cells(postnatal_cortex)[which(postnatal_cortex$replicates == "rep3rep4")]
slct_WT1_WT3 = sample(WT1_WT3, size = 5000)
slct_WT2_WT4 = sample(WT2_WT4, size = 5000)
downsampled_cells = subset(postnatal_cortex, cells = c(slct_WT1_WT3, slct_WT2_WT4))

DimPlot_scCustom(seurat_object = postnatal_cortex, split.by = 'orig.ident', pt.size = 0.8, label = FALSE)
ggplot2::ggsave(glue("{base_path}/P150_replicate.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)
############################################################
## Perform DEG analysis WT cells from the WT mouse cortex ##
############################################################
celltype_groups <- unique(downsampled_cells@meta.data$broad_class)
deg_results <- list()

for (celltype_group in celltype_groups) {
  cat("Performing DEG analysis for", celltype_group, "\n")
  
  # Subset cells based on celltype
  celltype_subset <- subset(downsampled_cells, subset = broad_class == celltype_group)
  
  # Get expression info
  expr <- as.matrix(GetAssayData(celltype_subset))
  
  # Filter out genes that are 0 for every cell
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad, ]
  
  logcpm <- cpm(expr, prior.count = 2, log = TRUE)
  mm <- model.matrix(~0 + replicates, data = celltype_subset@meta.data)
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Extract DEG results
  contrasts <- makeContrasts(c(replicatesrep1rep2) - c(replicatesrep3rep4), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contrasts)
  tmp <- eBayes(tmp)
  top_table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
  
  # Store DEG results for this celltype group
  deg_results[[celltype_group]] <- top_table
}
################################Test
for (celltype_group in celltype_groups) {
  cat("Performing DEG analysis for", celltype_group, "\n")
  
  # Subset cells based on celltype
  celltype_subset <- subset(downsampled_cells, subset = broad_class == celltype_group)
  
  # Get expression info
  expr <- as.matrix(GetAssayData(celltype_subset))
  
  # Filter out genes that are 0 for every cell
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad, ]
  
  # Check variation in gene expression
  if (nrow(expr) < 2) {
    cat("Not enough genes with variation for", celltype_group, "\n")
    next  # Skip to the next cell type
  }
  
  # Filter low-expressed genes
  expr <- expr[rowSums(expr) >= threshold, ]
  
  logcpm <- cpm(expr, prior.count = 2, log = TRUE)
  mm <- model.matrix(~0 + replicates, data = celltype_subset@meta.data)
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Extract DEG results
  contrasts <- makeContrasts(c(replicatesrep1rep2) - c(replicatesrep3rep4), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contrasts)
  tmp <- eBayes(tmp)
  top_table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
  
  # Store DEG results for this celltype group
  deg_results[[celltype_group]] <- top_table
}
################################Test End
deg_results$Glutamatergic$broad_class <- 'Glutamatergic'
deg_results$Glutamatergic$DEG_Test <- 'WT cells from P150 WT vs WT cells from P150 WT'
#MUTvsWT <- rbind(deg_results$P30, deg_results$P60, deg_results$P150)
#MUTvsWT$SYMBOL <- rownames(MUTvsWT)
#write.csv(MUTvsWT, file = glue('{base_path}/mutvswt_DEG_Glut.csv'))
####################
## Create volcano ##
####################
top.table <- deg_results$Glutamatergic
top.table$Gene <- rownames(top.table)
top.table$diffexpressed <- 'NO'
top.table$diffexpressed[top.table$logFC > 0 & top.table$adj.P.Val < 0.05] <- 'UP'
top.table$diffexpressed[top.table$logFC < 0 & top.table$adj.P.Val < 0.05] <- 'DOWN'
top.table$diffexpressed[top.table$adj.P.Val > 0.05] <- 'Not Sig'
top.table$delabel <- NA
thresh = head(arrange(top.table, adj.P.Val), 10)$adj.P.Val[10]
top.table$delabel[top.table$adj.P.Val <=thresh] <-(top.table$Gene[top.table$adj.P.Val<=thresh])
# Count the number of UP, DOWN, and Not Sig genes
count_genes <- table(top.table$diffexpressed)

# Print the counts
cat("Number of UP genes:", count_genes["UP"], "\n")
cat("Number of DOWN genes:", count_genes["DOWN"], "\n")
cat("Number of Not Sig genes:", count_genes["Not Sig"], "\n")
# Volcano Plot
# Count the number of UP, DOWN, and Not Sig genes
count_genes <- table(top.table$diffexpressed)

# Create the plot with the information added to the title
ggplot(data = top.table, aes(x = logFC, y = -log2(adj.P.Val), col = diffexpressed, label = delabel)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf, box.padding = 0.8) +
  scale_color_manual(values = c('blue','black','red')) +
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
    legend.text = element_text(size = 14, face = "bold")
  ) +
  labs(title = 'P150 replicates test',
       subtitle = paste('UP:', count_genes["UP"], '  ',
                        'DOWN:', count_genes["DOWN"], '  ',
                        'Not Sig:', count_genes["Not Sig"]))
ggplot2::ggsave(glue("{base_path}/Vol_gaba_reps_cort1_3.pdf"),
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
pdf(glue("{base_path}/venn_Glutamatergic.pdf"))
temp <- venn.diagram(
  x = list(
    P30 = sig_genes_P30,
    P60 = sig_genes_P60,
    P150 = sig_genes_P150
  ),
  category.names = c("P30", "P60", "P150"),
  main = 'Glutamatergic DEGs from MUT cells from MUT females and WT cells from WT females ',
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
    openxlsx::write.xlsx(file=glue::glue("{base_path}/GLUT_P60_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Mm.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{base_path}/GLUT_P60_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("{base_path}/GLUT_P60_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
  })
print(glue::glue("The pipeline has finished for samples"))
