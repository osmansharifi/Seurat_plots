library(topGO)
library(org.Mm.eg.db)
library(foreach)
library(glue)
library(ggplot2)
library(scales)
library(tidyverse)

# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

## Paths
#Limma_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_E18_WB/"
#Limma_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_P30_CORT/"
#Limma_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_P60_CORT/"
Limma_DEG_dir <- "~/GitHub/snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_P120_CORT/"

figure_path <- "~/GitHub/snRNA-seq-pipeline/figures/go_analysis/enrichment_scores/"


## Lists
cell_types <- list("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo") 
topgo_ontologies <- list("BP", "CC", "MF")


## Other variables
#metadata_info_concise <- "M_MUT_and_WT_M_E18_WB"
#metadata_info_concise <- "M_MUT_and_WT_M_P30_CORT"
#metadata_info_concise <- "M_MUT_and_WT_M_P60_CORT"
metadata_info_concise <- "M_MUT_and_WT_M_P120_CORT"

#metadata_info_expanded <- "Male, E18, Whole Brain"
#metadata_info_expanded <- "Male, P30, Cortex"
#metadata_info_expanded <- "Male, P60, Cortex"
metadata_info_expanded <- "Male, P120, Cortex"
################################################################################

typeof(geneList)
head(geneList)


for (cell_type in cell_types){
  # Read in significant DEGs per cell type identified by Limma
  signif_DEGs <- read.csv(file = glue(Limma_DEG_dir, cell_type, "_", metadata_info_concise, "_Limma_DEG.csv"))
  # Define geneList
  geneList <- as.vector(signif_DEGs[,c(1,6)])
  geneList <- geneList %>% remove_rownames %>% column_to_rownames(var="X")
  geneList$adj.P.Val <- as.numeric(geneList$adj.P.Val)
  
  
  
  for (ont in topgo_ontologies){
    # Create topGOdata object
    assign(glue('GOdata{ont}'), new("topGOdata",
                                    ontology = ont,
                                    allGenes = geneList,
                                    geneSelectionFun = function(x)(x == 1),
                                    annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol"))
  }
  godata_types <- list(GOdataBP, GOdataCC, GOdataMF)
  godata_names <- list("GOdataBP", "GOdataCC", "GOdataMF")
  foreach(GOdata = godata_types, godata_name = godata_names) %do% {
    # Test for enrichment using Fisher's Exact Test and visualize GO terms
    resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
    goEnrichment <- GenTable(
      GOdata,
      Fisher = resultFisher,
      orderBy = "Fisher",
      topNodes = 20,
      numChar = 50)
    goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
    # Filter terms for Fisher p<0.05
    goEnrichment <- goEnrichment[goEnrichment$Fisher < 0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","Fisher")]
    ntop <- 20
    ggdata <- goEnrichment[1:ntop,]
    ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
    ggplot(ggdata,
           aes(x = Term, y = -log10(Fisher), size = -log10(Fisher), fill = -log10(Fisher))) +
      expand_limits(y = 1) +
      geom_point(shape = 21) +
      scale_size(range = c(2.5,12.5)) +
      scale_fill_continuous(low = 'royalblue', high = 'red4') +
      xlab('') + ylab('Enrichment score') +
      labs(
        title = glue('GO Analysis using ', godata_name, ' for ', cell_type),
        subtitle = glue('Top 20 terms ordered by Fisher Exact p-value for {metadata_info_expanded}'),
        caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
      geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
                 linetype = c("dotted", "longdash", "solid"),
                 colour = c("black", "black", "black"),
                 size = c(0.5, 1.5, 3)) +
      theme_bw(base_size = 24) +
      theme(
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
        # Legend
        legend.key = element_blank(), # removes the border
        legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
        legend.text = element_text(size = 14, face = "bold"), # Text size
        title = element_text(size = 14, face = "bold")) +
      coord_flip()
    ggplot2::ggsave(glue('{figure_path}{cell_type}_{metadata_info_concise}_{godata_name}_Fisher.pdf'),
                    device = NULL,
                    height = 8.5,
                    width = 12)
  }
}

