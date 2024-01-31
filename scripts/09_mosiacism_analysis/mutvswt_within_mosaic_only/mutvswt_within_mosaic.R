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
##################
## Load samples ##
##################
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
load(glue('{base_path}/all.female.cortex.parsed.RData'))
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

all.female.cortex <- subset(x = all.female.cortex, subset = Condition == 'MUTANT')

# Perform DEG analysis between the Mut and WT cells within the mosaic brains
MUT_from_MUT = Cells(all.female.cortex)[which(all.female.cortex$Mecp2_allele == "Mecp2_MUT")]
WT_from_MUT = Cells(all.female.cortex)[which(all.female.cortex$Mecp2_allele == "Mecp2_WT")]
slct_MUT_from_MUT = sample(MUT_from_MUT, size = 539)
slct_WT_from_MUT = sample(WT_from_MUT, size = 539)
subset_cell_nonautonomous = subset(all.female.cortex, cells = c(slct_MUT_from_MUT, slct_WT_from_MUT))
subset_cell_nonautonomous <- subset(x = subset_cell_nonautonomous, subset = broad_class == 'GABAergic')
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
  mm <- model.matrix(~0 + Mecp2_allele, data = age_subset@meta.data)
  y <- voom(expr, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  
  # Extract DEG results
  contrasts <- makeContrasts(c(Mecp2_alleleMecp2_MUT) - c(Mecp2_alleleMecp2_WT), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contrasts)
  tmp <- eBayes(tmp)
  top_table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
  
  # Store DEG results for this age group
  deg_results[[age_group]] <- top_table
}

# Access DEG results for each age group
for (age_group in age_groups) {
  cat("DEG analysis results for", age_group, "\n")
  
  deg_table <- deg_results[[age_group]]
  num_degs <- length(which(deg_table$adj.P.Val < 0.05))
  
  print(num_degs)
  
  # Additional analysis if needed
  # summary(decideTests(tmp))
}

top.table <- deg_results$P30
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
  scale_color_manual(values = c('black', 'red')) +
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
  labs(title = 'P30 GABAergic MUT vs WT cells within mosaic brains')
ggplot2::ggsave(glue("{base_path}mutvswt_within_mosaic_only/GABA_P30_MUTvsWT_withinfemales.pdf"),
                device = NULL,
                height = 8.5,
                width = 12)

#Write csv files
deg_results$P30$Time_point <- "P30"
deg_results$P60$Time_point <- "P60"
deg_results$P150$Time_point <- "P150"
MUTvsWT <- rbind(deg_results$P30, deg_results$P60, deg_results$P150)
MUTvsWT$SYMBOL <- rownames(MUTvsWT)
write.csv(MUTvsWT, file = glue('{base_path}mutvswt_within_mosaic_only/MUTvsWT_within_mosaic_DEGs_GABA.csv'))

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
pdf(glue("{base_path}/broad_group_analysis/venn_Glutamatergic.pdf"))
temp <- venn.diagram(
  x = list(
    P30 = sig_genes_P30,
    P60 = sig_genes_P60,
    P150 = sig_genes_P150
  ),
  category.names = c("P30", "P60", "P150"),
  main = 'Glutamatergic DEGs from WT cells from WT females and WT cells from mosaic females ',
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
write.csv(as.data.frame(go.obj), file = glue('{base_path}/broad_group_analysis/geneoverlap_gaba_3timepoints.txt'))
DEGs = top.table
# Top DEGs

DEGs <- deg_results$P30 %>%
  tibble::rownames_to_column() %>%
  tibble::as_tibble() %>%
  dplyr::rename(SYMBOL = rowname) %>%
  dplyr::mutate(FC = dplyr::case_when(logFC >0 ~ 2^logFC,
                                      logFC <0 ~ -1/(2^logFC))) %>%
  dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %T>%
  openxlsx::write.xlsx(file=glue::glue("{base_path}mutvswt_within_mosaic_only/DEGs.xlsx")) %>%
  dplyr::filter(P.Value < 0.05) %T>%
  openxlsx::write.xlsx(file=glue::glue("{base_path}mutvswt_within_mosaic_only/sig_DEGs.xlsx"))

enrichR:::.onAttach()
source(glue('{base_path}/GO_ploting_functions.R'))
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
    openxlsx::write.xlsx(file=glue::glue("{base_path}mutvswt_within_mosaic_only/GABA_P30_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Mm.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{base_path}mutvswt_within_mosaic_only/GABA_P30_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("{base_path}mutvswt_within_mosaic_only/GABA_P30_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
  })
