## Filtering genes into highly and lowly expressed for downstream analyses

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "limma", "edgeR", "ggplot2", "ggpubr", "viridis", "recommenderlab")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

setwd("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/Mecp2e1_parsed")

load("DEanalysis_Mecp2_parsing_00.RData")

## exploring data

sum = design %>%
  dplyr::group_by(sex, age, genotype, Mecp2e1_expression) %>%
  dplyr::summarize(n())

WTE18M = design %>%
  dplyr::filter(sex == "Male", age=="E18", genotype=="WT", Mecp2e1_expression=="MUTANT")


## removing E18 mice for age combined analyses

expr_matrix_noE18 <- expr_matrix[,which(design$age != "E18")]

design <- design %>%
  dplyr::filter(age != "E18")

design$age = droplevels(design$age)

## creating count matrices split by cell type

design$cell_type <- droplevels(design$cell_type)

cell_types = as.factor(levels(design$cell_type))

expr_matrix_list = lapply(levels(cell_types), function(x) {
  expr_matrix[,design$cell_ID[which(design$cell_type==x)]]
})

names(expr_matrix_list) = as.character(cell_types)

## Creating DGEList

DGEList = lapply(cell_types, function(x) {
  DGEList(expr_matrix_list[[x]])
})

names(DGEList) = cell_types

### Filtering DGEList by highly and lowly expressed genes. Highly expressed = at least 1 CPM in more than 25% of cells

DGEListCPM = lapply(cell_types, function(x){
  cpm(DGEList[[x]])
})

names(DGEListCPM) = cell_types

highly_expr_genes = lapply(cell_types, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})

names(highly_expr_genes) = cell_types

data = sapply(cell_types, function(x) {length(highly_expr_genes[[x]])})
data.frame1 = data.frame(cell_type = cell_types, num_highly_expressed_genes=data, cut_off = rep("25%", times=length(cell_types)))

highly_expr_genes = lapply(cell_types, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.50*ncol(DGEList[[x]]$counts))=="TRUE")
})

names(highly_expr_genes) = cell_types

data = sapply(cell_types, function(x) {length(highly_expr_genes[[x]])})
data.frame2 = data.frame(cell_type = cell_types, num_highly_expressed_genes=data, cut_off = rep("50%", times=length(cell_types)))

data.frame.comb = rbind(data.frame1, data.frame2)

pdf(file="Gene_Filtering_cutoffs_25_50.pdf", height=8.5, width=11)
ggplot(data=data.frame.comb, aes(x=cell_type, y=num_highly_expressed_genes)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  facet_grid(~cut_off) +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="NULL") +
  ggtitle("Genes with at least 1 CPM in more than either 25% or 50% of cells") 
dev.off()

# Using 25% cutoff
highly_expr_genes = lapply(cell_types, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})
names(highly_expr_genes) = cell_types

lowly_expr_genes = lapply(cell_types, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) < 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})
names(lowly_expr_genes) = cell_types

DGEList_high = lapply(cell_types, function(x){
  DGEList[[x]][highly_expr_genes[[x]],] %>%
    calcNormFactors()
})

DGEList_low = lapply(cell_types, function(x){
  DGEList[[x]][lowly_expr_genes[[x]],] %>%
    calcNormFactors()
})

names(DGEList_high) = cell_types
names(DGEList_low) = cell_types

save(DGEList, DGEList_high, DGEList_low, design, cell_types, file="DEanalysis_Mecp2_parsing_age_combined_01.RData")



