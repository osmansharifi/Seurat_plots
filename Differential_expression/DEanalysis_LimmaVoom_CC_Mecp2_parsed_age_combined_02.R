### Differential Expression of wt-expressing cells by Mecp2-e1 genotype with all ages combined except E18
### LimmaVoom CC method

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Glimma", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis", "RColorBrewer", "BiocParallel", "DEsingle", "enrichR", "DMRichR", "org.Mm.eg.db", "AnnotationDbi")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))
enrichR:::.onAttach()

setwd("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/Mecp2e1_parsed")

load(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/Mecp2e1_parsed/DEanalysis_Mecp2_parsing_age_combined_01.RData.RData"))

experiment = "genotype_age_combined" #change to reflect age x genotype analysis directory

dir.create("limmaVoomCC")
setwd("limmaVoomCC")

dir.create(glue::glue("{experiment}"))
setwd(glue::glue("{experiment}"))

cell_types_all = names(DGEList_high)

### Looking at only wild-type expressing cells

DGEList_high_new <- vector(mode='list', length=length(cell_types_all))
names(DGEList_high_new) = cell_types_all
  
for(i in cell_types_all) {
  design_sub = design %>%
    dplyr::filter(cell_type == i)
  DGEList_high_new[[i]] <- DGEList_high[[i]][,which(design_sub$Mecp2e1_expression == "WT")]
}

for (i in cell_types_all) {
  DGEList_high_new[[i]]$samples$group = ifelse(grepl("WT", colnames(DGEList_high_new[[i]]$counts))=="TRUE", "WT", "MUTANT")
}

cell_types_all = cell_types_all[-c(2)] # removing Endo bc no cells

for (i in cell_types_all) {
  
  dir.create(glue::glue("{i}"))
  dir.create(glue::glue("{i}/QC"))
  dir.create(glue::glue("{i}/plotData"))
  
  print(glue::glue("Normalizing {i} cells"))
  
  #subsetting design matrix by cell type and filtering out mutant-expressing cells
  
  design.new <- design %>%
    dplyr::filter(cell_type == i) %>%
    dplyr::filter(Mecp2e1_expression == "WT")
  
  design.new$genotype = factor(design.new$genotype, levels=c("WT", "MUTANT"))
  
  # Raw density of log-CPM values
  
  logCPM <- cpm(DGEList_high_new[[i]], log=TRUE)
  nsamples <- ncol(DGEList_high_new[[i]])
  col <- brewer.pal(nsamples, "Paired")
  
  pdf(glue::glue("{i}/density_plot.pdf"), height=8.5, width=5.5)
  
  plot(density(logCPM[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="Density Plot of Highly Expressed Genes", xlab="Log-cpm")
  for(j in 2:nsamples){
    den <- density(logCPM[,j])
    lines(den$x, den$y, col=col[j], lwd=2)
  }
  legend("topright", colnames(DGEList_high_new[[i]]$counts), text.col=col, bty="n", cex=0.5)
  
  dev.off()
  
  Glimma::glMDSPlot(DGEList_high_new[[i]],
                    groups = design.new$genotype,
                    path = getwd(),
                    folder = "interactivePlots",
                    html = glue::glue("{i}_MDS-Plot"),
                    launch = FALSE)
  
  
  mm <- model.matrix(~genotype + age,
                     data = design.new) #unable to correct for %mito and cell cycle due to insufficient numbers of cells that span all of these conditions
  
  
  # Voom
  pdf(glue::glue("{i}/QC/voom_mean-variance_trend.pdf"), height=8.5, width=11)
  voomLogCPM <- voom(DGEList_high_new[[i]],
                     mm,
                     plot=TRUE)
  dev.off()
  
  correlations <- duplicateCorrelation(voomLogCPM,
                                       mm,
                                       block = design.new$sample_ID)
  
  correlations <- correlations$consensus.correlation
  
  # Boxplots of logCPM values before and after voom normalization
  pdf(glue::glue("{i}/QC/normalization_boxplots.pdf"), height=8.5, width=11)
  par(mfrow=c(1,2))
  
  boxplot(logCPM, las=2, col=col, main="")
  title(main="A. Unnormalized data", ylab="Log-cpm")
  
  boxplot(voomLogCPM$E, las=2, col=col, main="")
  title(main="B. Normalized data", ylab="Log-cpm")
  
  dev.off()
  
  # Fitting linear models in limma
  
  print(glue::glue("Testing {i} cells for differential expression"))
  
  fit <- lmFit(voomLogCPM,
               mm,
               correlation = correlations,
               block = design.new$sample_ID)
  
  head(coef(fit))
  
  # Save normalized expression values
  voomLogCPM$E %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "geneSymbol") %>%
    openxlsx::write.xlsx(glue::glue("{i}/plotData/voomLogCPM.xlsx"))
  
  # Create DEG List 
  
  print(glue::glue("Creating DEG list of {i} cells for genotype"))
  
  efit <- fit %>%
    contrasts.fit(coef="genotypeMUTANT") %>%
    eBayes()
  
  pdf(glue::glue("{i}/QC/final_model_mean-variance_trend.pdf"), height=8.5, width=11)
  
  plotSA(efit, main=glue::glue("Final model {i}: Mean-variance trend"))
  
  dev.off()
  
  Glimma::glimmaMA(efit,
                   dge = DGEList_high_new[[i]],
                   path = getwd(),
                   html = glue::glue("interactivePlots/{i}_MDA-Plot.html"),
                   launch = FALSE)
  
  # Top DEGs
  
  DEGs <- efit %>%
    topTable(sort.by = "P", n = Inf) %>%
    tibble::rownames_to_column() %>%
    tibble::as_tibble() %>%
    dplyr::rename(SYMBOL = rowname) %>%
    dplyr::mutate(FC = dplyr::case_when(logFC >0 ~ 2^logFC,
                                        logFC <0 ~ -1/(2^logFC))) %>%
    dplyr::select(SYMBOL, FC, logFC, P.Value, adj.P.Val, AveExpr, t, B) %T>%
    openxlsx::write.xlsx(file=glue::glue("{i}/DEGs.xlsx")) %>%
    dplyr::filter(P.Value < 0.05) %T>%
    openxlsx::write.xlsx(file=glue::glue("{i}/sig_DEGs.xlsx"))
  
  print(glue::glue("GO and Pathway analysis of {i} cells"))
  
  enrichR:::.onAttach()
  
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
  print(glue::glue("The pipeline has finished for {i} cells"))

  
}







