packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Glimma", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis", "RColorBrewer", "BiocParallel", "DEsingle", "enrichR", "DMRichR", "org.Hs.eg.db", "AnnotationDbi")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

enrichR:::.onAttach()

param <- SnowParam(2, "SOCK", progressbar = TRUE)
register(param)

s.obj.name = "human_rett_cort_filt" #change this to the name of the Seurat object you started with

setwd(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/{s.obj.name}"))

load("DEanalysis_01.RData")
load("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rettcort_labeled/DEanalysis_01.RData")
dir.create("limmaVoomCC")
setwd("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/limmaVoomCC")

polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")
#DGEList_high = DGEList_high[-c(which(names(DGEList_high)=="Non-neuronal"))] #removing non-neuronal cells due to lack of sufficient cell numbers across genotypes

levels(design$cell_type) <- c("Astro", "Endo", "L2_3_IT", "L5", "L5", "L5_6", "L6", "L6", "L6", "L6", "Lamp5","Micro","Oligo", "OPC", "Pvalb", "Sncg", "Sst", "Sst", "Vip", "Micro")
DGEList_high$L2_3_IT <- DGEList_high$`L2/3 IT`

cell_types_all = names(DGEList_high)

#cell_types_all = cell_types_all[-(c(which(grepl("Sncg-activated", cell_types_all)=="TRUE")))] #removing Sncg-activated from analysis bc not enough cells

cell_types_all = cell_types_all[c(3:3)] # not enough cells for these cell types

for (i in cell_types_all) {
  DGEList_high[[i]]$samples$group = ifelse(grepl("CTRL", colnames(DGEList_high[[i]]$counts))=="TRUE", "CTRL", "RTT")
}

cell_types_all = c("Astro", "Endo", "L2_3_IT", "L5", "L5", "L5_6", "L6", "L6", "L6", "L6", "Lamp5","Micro","Oligo", "OPC", "Pvalb", "Sncg", "Sst", "Sst", "Vip", "Micro")

cell_types_all = cell_types_all[c(2:18)]

for (i in cell_types_all) {
  
  dir.create(glue::glue("{i}"))
  dir.create(glue::glue("{i}/QC"))
  dir.create(glue::glue("{i}/plotData"))
  
  print(glue::glue("Normalizing {i} cells"))
  
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
  
  design.new$genotype = factor(design.new$genotype, levels=c("CTRL", "RTT"))
  
  # Raw density of log-CPM values
  
  L <- mean(DGEList_high[[i]]$samples$lib.size)*1e-6
  M <- median(DGEList_high[[i]]$samples$lib.size)*1e-6
  
  logCPM <- cpm(DGEList_high[[i]], log=TRUE)
  nsamples <- ncol(DGEList_high[[i]])
  #col <- brewer.pal(nsamples, "Paired")
  col <- polychrome_palette
  
  pdf(glue::glue("{i}/density_plot.pdf"), height=8.5, width=5.5)
  
  plot(density(logCPM[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="Density Plot of Highly Expressed Genes", xlab="Log-cpm")
  for(j in 2:nsamples){
    den <- density(logCPM[,j])
    lines(den$x, den$y, col=col[j], lwd=2)
  }
  legend("topright", colnames(DGEList_high[[i]]$counts), text.col=col, bty="n", cex=0.5)
  
  dev.off()
  
  #Glimma::glMDSPlot(DGEList_high[[i]],
                    #groups = design.new$genotype,
                   # path = getwd(),
                    #folder = "interactivePlots",
                    #html = glue::glue("{i}_MDS-Plot"),
                   # launch = FALSE)
  
  
  mm <- model.matrix(~genotype + cell_cycle + percent.mito,
                     data = design.new)
  
  
  # Voom
  pdf(glue::glue("{i}/QC/voom_mean-variance_trend.pdf"), height=8.5, width=11)
  voomLogCPM <- voom(DGEList_high[[i]],
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
    contrasts.fit(coef="genotypeRTT") %>%
    eBayes()
  
  pdf(glue::glue("{i}/QC/final_model_mean-variance_trend.pdf"), height=8.5, width=11)
  
  plotSA(efit, main=glue::glue("Final model {i}: Mean-variance trend"))
  
  dev.off()
  
  #Glimma::glimmaMA(efit,
                  # dge = DGEList_high[[i]],
                  # path = getwd(),
                   #html = glue::glue("interactivePlots/{i}_MDA-Plot.html"),
                   #launch = FALSE)
  
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
                      annoDb = "org.Hs.eg.db",
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
