#### Differential Expression Analysis of scRNA-seq Data - 02 - limma voom pseudobulk (PB) analysis ####

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Glimma", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis", "RColorBrewer", "BiocParallel", "DEsingle", "enrichR", "DMRichR", "org.Mm.eg.db", "AnnotationDbi")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

enrichR:::.onAttach()

param <- SnowParam(2, "SOCK", progressbar = TRUE)
register(param)

s.obj.name = "rett_E18_with_labels_proportions" #change this to the name of the Seurat object you started with

setwd(glue::glue("/Users/karineier/Documents/scRNA-seq-differential-expression/{s.obj.name}"))

load("DEanalysis_01.RData")

dir.create("limmaVoomPB")
setwd("limmaVoomPB")

DGEList = DGEList[c(19:30)] # no subsetting by activated vs. unactivated due to lack of sufficient cells across samples

cell_types_all = names(DGEList)

# Aggregating gene counts across cells within each cell type

DGEList_aggregate = lapply(cell_types_all, function(i){
  
  design.new <- design %>%
    dplyr::filter(cell_type == i)
  
  design.new$sample_ID = droplevels(design.new$sample_ID)
  
  new.mat = sapply(1:length(levels(design.new$sample_ID)), function(x){
    rowSums(DGEList[[i]]$counts[,which(design.new$sample_ID==as.character(levels(design.new$sample_ID)[x]))])
  })
  colnames(new.mat) = levels(design.new$sample_ID)
  
  return(new.mat)
})

names(DGEList_aggregate) = cell_types_all

for (i in cell_types_all) {
  
  dir.create(glue::glue("{i}"))
  dir.create(glue::glue("{i}/QC"))
  dir.create(glue::glue("{i}/plotData"))
  
  print(glue::glue("Normalizing {i} cells"))
  
  design.agg <- data.frame(genotype=c(ifelse(grepl("WT", colnames(DGEList_aggregate[[i]]))=="TRUE", "WT", "MUTANT")))
  
  design.agg$genotype = factor(design.agg$genotype, levels=c("WT", "MUTANT"))
  
  #create DGEList
  DGEList_aggregate[[i]] = DGEList_aggregate[[i]] %>%
    DGEList() %>%
    calcNormFactors()
  
  # Raw density of log-CPM values
  
  L <- mean(DGEList_aggregate[[i]]$samples$lib.size)*1e-6
  M <- median(DGEList_aggregate[[i]]$samples$lib.size)*1e-6
  
  logCPM <- cpm(DGEList_aggregate[[i]], log=TRUE)
  logCPM.cutoff <- log2(10/M + 2/L)
  nsamples <- ncol(DGEList_aggregate[[i]])
  col <- brewer.pal(nsamples, "Paired")
  
  pdf(glue::glue("{i}/density_plot.pdf"), height=8.5, width=11)
  par(mfrow=c(1,2))
  
  plot(density(logCPM[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  abline(v=logCPM.cutoff, lty=3)
  for(j in 2:nsamples){
    den <- density(logCPM[,j])
    lines(den$x, den$y, col=col[j], lwd=2)
  }
  legend("topright", colnames(DGEList_aggregate[[i]]$counts), text.col=col, bty="n", cex=0.5)
  
  # Filter genes with low expression
  
  rawCount <- dim(DGEList_aggregate[[i]])
  
  keep.exprs <- filterByExpr(DGEList_aggregate[[i]],
                             group = design.agg$genotype,
                             lib.size = DGEList_aggregate[[i]]$samples$lib.size)
  
  DGEList_aggregate[[i]] <- DGEList_aggregate[[i]][keep.exprs,, keep.lib.sizes=FALSE] %>%
    calcNormFactors()
  
  filterCount <- dim(DGEList_aggregate[[i]])
  
  print(glue::glue("{100 - round((filterCount[1]/rawCount[1])*100)}% of genes were filtered from {rawCount[2]} samples, \\
                   where there were {rawCount[1]} genes before filtering and {filterCount[1]} genes after filtering for {i} cells"))
  
  logCPM <- cpm(DGEList_aggregate[[i]], log=TRUE)
  plot(density(logCPM[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="B. Filtered data", xlab = "Log-cpm")
  abline(v=logCPM.cutoff, lty=3)
  for(j in 2:nsamples){
    den<- density(logCPM[,j])
    lines(den$x, den$y, col=col[j], lwd=2)
  }
  legend("topright", colnames(DGEList_aggregate[[i]]$counts), text.col=col, bty="n", cex=0.5)
  dev.off()
  
  Glimma::glMDSPlot(DGEList_aggregate[[i]],
                    groups = design.agg$genotype,
                    path = getwd(),
                    folder = "interactivePlots",
                    html = glue::glue("{i}_MDS-Plot"),
                    launch = FALSE)
  
  mm <- model.matrix(~genotype,
                     data = design.agg)
  
  # Voom
  pdf(glue::glue("{i}/QC/voom_mean-variance_trend.pdf"), height=8.5, width=11)
  voomLogCPM <- voom(DGEList_aggregate[[i]],
                     mm,
                     plot=TRUE)
  dev.off()
  
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
               mm)
  
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


