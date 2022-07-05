### Differential Expression Analysis of scRNA-seq Data - 00 - Splitting Seurat Object ###

# Splitting Seurat Object by Activated and Unactivated Neurons

# set up 

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "Matrix.utils", "edgeR", "ggplot2")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

s.obj.names = c("all_female_P30", "all_female_P60", "all_female_P150") #change this to the name of the Seurat object you're working with

for (i in s.obj.names) {
  load(glue::glue("/Users/karineier/Documents/Mecp2/scRNA-seq/{i}.RData"))
}

# merge count matrices keeping all genes present in at least one Seurat object with the all.x = TRUE and all.y = TRUE arguments

count.matrix.1 = merge.Matrix(all_female_P30@assays$RNA@counts, all_female_P60@assays$RNA@counts, 
                              by.x = rownames(all_female_P30@assays$RNA@counts), by.y = rownames(all_female_P60@assays$RNA@counts), 
                              all.x=TRUE, all.y=TRUE)
count.matrix.2 = merge.Matrix(count.matrix.1, all_female_P150@assays$RNA@counts, 
                              by.x = rownames(count.matrix.1), by.y = rownames(all_female_P150@assays$RNA@counts), 
                              all.x=TRUE, all.y=TRUE)

dir.create(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/females_by_age"))

setwd(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/females_by_age"))

# making design data frame

design = data.frame(cell_type = c(all_female_P30$predicted.id, all_female_P60$predicted.id, all_female_P150$predicted.id),
                    sample_ID = c(all_female_P30$orig.ident, all_female_P60$orig.ident, all_female_P150$orig.ident),
                    cell_cycle = c(all_female_P30$Phase, all_female_P60$Phase, all_female_P150$Phase),
                    cell_ID = c(colnames(all_female_P30@assays$RNA@counts), colnames(all_female_P60@assays$RNA@counts), colnames(all_female_P150@assays$RNA@counts)),
                    genotype = c(ifelse(grepl("WT", all_female_P30$orig.ident)=="TRUE", "WT", "MUTANT"),
                                 ifelse(grepl("WT", all_female_P60$orig.ident)=="TRUE", "WT", "MUTANT"),
                                 ifelse(grepl("WT", all_female_P150$orig.ident)=="TRUE", "WT", "MUTANT")),
                    percent.mito = c(all_female_P30$percent.mito, all_female_P60$percent.mito, all_female_P150$percent.mito),
                    age = c(rep(30, times=ncol(all_female_P30@assays$RNA@counts)), 
                            rep(60, times=ncol(all_female_P60@assays$RNA@counts)),
                            rep(150, times=ncol(all_female_P150@assays$RNA@counts)))
                    )

design = design %>%
  dplyr::mutate_if(is.character, as.factor)

design$genotype = factor(design$genotype, levels=c("WT", "MUTANT"))

# checking to make sure colnames of design and count matrix are the same - output should be TRUE

identical(as.character(design$cell_ID), as.character(colnames(count.matrix.2)))

## creating count matrices split by cell type

cell_types = as.factor(levels(design$cell_type))

expr_matrix_list = lapply(levels(cell_types), function(x) {
  count.matrix.2[,design$cell_ID[which(design$cell_type==x)]]
})

names(expr_matrix_list) = as.character(cell_types)

# create DGEList object

cell_types_all = names(expr_matrix_list)

DGEList = lapply(cell_types_all, function(x) {
  DGEList(expr_matrix_list[[x]])
})

names(DGEList) = cell_types_all

### Filtering DGEList by highly and lowly expressed genes. Highly expressed = at least 1 CPM in more than 25% of cells

DGEListCPM = lapply(cell_types_all, function(x){
  cpm(DGEList[[x]])
})

names(DGEListCPM) = cell_types_all

highly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})

names(highly_expr_genes) = cell_types_all

data = sapply(cell_types_all, function(x) {length(highly_expr_genes[[x]])})
data.frame1 = data.frame(cell_type = cell_types_all, num_highly_expressed_genes=data, cut_off = rep("25%", times=length(cell_types_all)))

highly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.50*ncol(DGEList[[x]]$counts))=="TRUE")
})

names(highly_expr_genes) = cell_types_all

data = sapply(cell_types_all, function(x) {length(highly_expr_genes[[x]])})
data.frame2 = data.frame(cell_type = cell_types_all, num_highly_expressed_genes=data, cut_off = rep("50%", times=length(cell_types_all)))

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
highly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) > 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})
names(highly_expr_genes) = cell_types_all

lowly_expr_genes = lapply(cell_types_all, function(x){
  which((rowSums(DGEListCPM[[x]]>=1, na.rm=T) < 0.25*ncol(DGEList[[x]]$counts))=="TRUE")
})
names(lowly_expr_genes) = cell_types_all

DGEList_high = lapply(cell_types_all, function(x){
  DGEList[[x]][highly_expr_genes[[x]],] %>%
    calcNormFactors()
})

DGEList_low = lapply(cell_types_all, function(x){
  DGEList[[x]][lowly_expr_genes[[x]],] %>%
    calcNormFactors()
})

names(DGEList_high) = cell_types_all
names(DGEList_low) = cell_types_all

save(DGEList, DGEList_high, DGEList_low, design, cell_types_all, file="DEanalysis_00-01_combined_by_time_filtered.RData")
