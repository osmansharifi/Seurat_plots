
#### Plotting log2(fold-change) of genes vs. gene length ####

packages <- c("edgeR", "tidyverse", "RColorBrewer", "org.Mm.eg.db", "AnnotationDbi", "EnhancedVolcano", "enrichR", "openxlsx", "glue", "Glimma", "DMRichR", "magrittr", "variancePartition", "UpSetR", "ComplexUpset", "EDASeq", "broom")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

s.obj.name = "all_female_P150" # replace with the name of the original Seurat object
s.obj.name = "mousecortexsca"
#setwd(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/{s.obj.name}"))
setwd(glue::glue("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/females/M_MUT_and_WT_F_P30_CORT"))
load("/Users/osman/Desktop/LaSalle_lab/Rett_Data/E18/total_e18_labeled.RData")

### Reading in DEsingle files

cell_types = list.dirs(glue::glue("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/females/M_MUT_and_WT_F_P30_CORT/limmaVoomCC"))
#cell_types = cell_types[-c(1, which(grepl("QC", cell_types)=="TRUE"), which(grepl("plotData", cell_types)=="TRUE"))]
#cell_types = cell_types[-c(which(grepl("", cell_types)=="TRUE"), 23)]
cell_types = cell_types[-1]
cell_types

cell_types.all = cell_types[-c(which(grepl("-activated", cell_types) == "TRUE"), which(grepl("-not", cell_types)=="TRUE"))]

list <- lapply(cell_types.all, function(cellType) {
  openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
})

cell_types = gsub(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/{s.obj.name}/DEsingle/"), "", cell_types.all)

names(list) = cell_types

### Reading in limmaCC files

cell_types.l = list.dirs(glue::glue("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/females/M_MUT_and_WT_F_P30_CORT/limmaVoomCC"))
cell_types.l = cell_types.l[-c(1, which(grepl("QC", cell_types.l)=="TRUE"), which(grepl("plotData", cell_types.l)=="TRUE"), which(grepl("interactivePlots", cell_types.l)=="TRUE"))]
cell_types.l

cell_types.l = cell_types.l[-c(which(grepl("-activated", cell_types.l) == "TRUE"), which(grepl("-not", cell_types.l)=="TRUE"))]

list.l <- lapply(cell_types.l, function(cellType) {
  openxlsx::read.xlsx(glue::glue("{cellType}/DEGs.xlsx"), rowNames = TRUE)
})


cell_types.l = gsub(glue::glue("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/females/M_MUT_and_WT_F_P30_CORT/limmaVoomCC"), "", cell_types.l)

cell_types.l = gsub("/", "", cell_types.l)
names(list.l) = cell_types.l

list.l.select <- lapply(cell_types.l, function(x) {
  
  dplyr::select(list.l[[x]], logFC, P.Value) %>%
    dplyr::rename(log2FC = logFC) %>%
    dplyr::rename(pvalue = P.Value)
  
  })
names(list.l.select) = cell_types.l
### Get gene length and GC content for each gene in our DEG data

list.new <- vector(mode="list", length=length(cell_types.l))
names(list.new) = cell_types.l

for(i in cell_types.l) {
  list.l.select[[i]]$ENSEMBL <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                             keys=rownames(list.l.select[[i]]),
                                             keytype="SYMBOL",
                                             column="ENSEMBL")
  list.new[[i]] = list.l.select[[i]][!is.na(list.l.select[[i]]$ENSEMBL),]
  mat <- getGeneLengthAndGCContent(list.new[[i]]$ENSEMBL, "mm10", mode='org.db')
  list.new[[i]]$gene_length <- mat[,1]
  list.new[[i]]$GC_content <- mat[,2]
}


##### Is gene length associated with fold-change in mutant vs. wild-type mice? 

# Get log2(fold-change) as a measure of magnitude of change and subset data by up- and down-regulated genes

for(i in cell_types) {
  
  list.new[[i]]$log2FC = log2(list.new[[i]]$norm_foldChange)
  
}

# subset into up and down-regulated genes

list.up <- lapply(cell_types.l, function(x) {
  subset(list.new[[x]], log2FC>0)
})
names(list.up) = cell_types.l


list.down <- lapply(cell_types.l, function(x) {
  subset(list.new[[x]], log2FC<0)
})
names(list.down) = cell_types.l

## Associations between gene length and log2FC

# in up
assoc.list.up <- lapply(cell_types.l, function(x) {
  data = list.up[[x]] %>% 
    dplyr::filter(pvalue<0.05) # only looking at genes that are associated with Mecp2 genotype at p<0.05
  lm(log2FC ~ log2(gene_length), data=data)
})
names(assoc.list.up) = cell_types.l

assoc.df.up <- plyr::ldply(assoc.list.up, tidy, .id="Cell_Type")
assoc.df.up <- subset(assoc.df.up, !(term %in% "(Intercept)"))

write.xlsx(assoc.df.up, "Associations_btw_Gene_Length_and_FoldChange_Upregulated_genes.xlsx")

# in down
assoc.list.down <- lapply(cell_types.l, function(x) {
  list.down[[x]] <- list.down[[x]][which(is.infinite(list.down[[x]]$log2FC)==FALSE),] #this line removes genes where log2FC has infinite values
  data = list.down[[x]] %>% 
    dplyr::filter(pvalue<0.05) # only looking at genes that are associated with Mecp2 genotype at p<0.05
  #if(nrow(data) <1, print("NULL")) else (lm(log2FC ~ log2(gene_length), data =data)
  if (nrow(data) < 3 ){
    print("NULL")
  } else {
    lm(log2FC ~ log2(gene_length), data = data)
  }
})
names(assoc.list.down) = cell_types.l

assoc.df.down <- plyr::ldply(assoc.list.down, tidy, .id="Cell_Type")
assoc.df.down <- subset(assoc.df.down, !(term %in% "(Intercept)"))

write.xlsx(assoc.df.down, "Associations_btw_Gene_Length_and_FoldChange_Downregulated_genes.xlsx")

#### Plots

#list.up$Lamp5$Type <- as.character(list.up$Lamp5$Type)
#list.up$Lamp5$State <- as.character(list.up$Lamp5$State)

#list.up$Sncg$Type <- as.character(list.up$Sncg$Type)
#list.up$Sncg$State <- as.character(list.up$Sncg$State)

plotdata <- dplyr::bind_rows(list.up, .id="Cell_Type") %>%
  dplyr::filter(pvalue<0.1)

polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#66B0FFFF", "#3B00FBFF")
plotdata$Cell_Type = factor(plotdata$Cell_Type, levels = c("L2_3_IT", "L4", "L5", "L6", "Pvalb", "Vip", "Sst", "Sncg", "Lamp5", "Peri", "Endo", "Oligo", "Astro", "Non-neuronal"))

tiff(filename="gene_length_upregulated.tiff", res=300, height=5, width=6.5, units="in")
ggplot(data=plotdata, aes(x=log2(gene_length), y=log2FC, color=Cell_Type)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values=polychrome_palette) +
  theme_bw()
dev.off()

plotdata <- dplyr::bind_rows(list.down, .id="Cell_Type") %>%
  dplyr::filter(pvalue<0.1)

polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#66B0FFFF", "#3B00FBFF")
plotdata$Cell_Type = factor(plotdata$Cell_Type, levels = c("L2_3_IT", "L4", "L5", "L6", "Pvalb", "Vip", "Sst", "Sncg", "Lamp5", "Peri", "Endo", "Oligo", "Astro", "Non-neuronal"))

tiff(filename="gene_length_downregulated.tiff", res=300, height=5, width=6.5, units="in")
ggplot(data=plotdata, aes(x=log2(gene_length), y=log2FC, color=Cell_Type)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(values=polychrome_palette) +
  theme_bw()
dev.off()

save(list = ls(), file = "gene_length.rds")
