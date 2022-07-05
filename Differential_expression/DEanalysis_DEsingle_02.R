#### Differential Expression Analysis of scRNA-seq Data - 02 - DEsingle for lowly expressed genes ####

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis", "BiocParallel", "DEsingle", "enrichR", "DMRichR", "org.Mm.eg.db", "AnnotationDbi")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

enrichR:::.onAttach()

param <- SnowParam(2, "SOCK", progressbar = TRUE)
register(param)

s.obj.name = "rett_E18_with_labels_proportions" #change this to the name of the Seurat object you started with

setwd(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/{s.obj.name}"))

load("DEanalysis_01.RData")

dir.create("DEsingle")
setwd("DEsingle")

cell_types_all = cell_types_all[-c(1:18, 20, 26, 28)]

for (i in cell_types_all) {
  
  dir.create(glue::glue("{i}"))
  
  #subsetting design matrix by cell type
  activation = ifelse(grepl("activated", i)=="TRUE", "activated", ifelse(grepl("not", i)=="TRUE", "unactivated", "both"))
  
  cellType = sub("-.*", "", i)
  
  design.new <- design %>%
    dplyr::filter(cell_type == cellType)
  
  design.new <- if (activation == "activated") {
    dplyr::filter(design.new, activation.status == "activated")
  } else {
    if (activation == "unactivated") {
      dplyr::filter(design.new, activation.status == "unactivated")
    } else {
      dplyr::filter(design.new, activation.status == "activated" | activation.status == "unactivated")
    }
  }
  
  #checking to make sure columns of the DGEList counts matrix matches up with the order of the design matrix so that we appropriately model genotype
  stopifnot(colnames(DGEList_low[[i]]$counts)==as.character(design.new$cell_ID))
  
  
  ## DESingle 
  
  print(glue::glue("DEsingle for {i} cells"))
  
  results <- DEsingle(counts=DGEList_low[[i]]$counts, group=design.new$genotype)
  
  results.classified <- DEtype(results=results, threshold=0.1)
  
  pdf(file=glue::glue("{i}/histogram_pvalues.pdf"), height=8.5, width=11)
  hist(results.classified$pvalue, breaks=20)
  dev.off()
  
  results.final <- results.classified %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::rename(SYMBOL = rowname) %>%
    dplyr::select(SYMBOL, foldChange, norm_foldChange, pvalue, pvalue.adj.FDR, Type, State) %>%
    openxlsx::write.xlsx(file=glue::glue("{i}/DEGs.xlsx"))
  
  print(glue::glue("GO and Pathway analysis of {i} cells"))
  
  enrichR:::.onAttach()
  
  tryCatch({
    results.classified %>% 
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::rename(SYMBOL = rowname) %>%
      dplyr::filter(pvalue<0.05) %>%
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
      openxlsx::write.xlsx(file=glue::glue("{i}/enrichr.xlsx")) %>%
      DMRichR::slimGO(tool = "enrichR",
                      annoDb = "org.Mm.eg.db",
                      plots = FALSE) %T>%
      openxlsx::write.xlsx(file = glue::glue("{i}/rrvgo_enrichr.xlsx")) %>%
      DMRichR::GOplot() %>%
      ggplot2::ggsave(glue::glue("{i}/enrichr_plot.pdf"),
                      plot = .,
                      device = NULL,
                      height = 8.5,
                      width = 10) },
    error = function(error_condition) {
      print(glue::glue("ERROR: Gene Ontology pipe did not finish for {i} cells"))
    })
    
}


