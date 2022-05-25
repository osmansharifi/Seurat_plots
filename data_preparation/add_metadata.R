library(Seurat)
library(dplyr)
library(scCustomize)
library(patchwork)

#Load data
load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all.cortex.combined.RData")
counts_table <- read.table("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/data_preparation/master_sequence_alleler.txt", sep="\t", header=FALSE)
#headers for the txt file
names(counts_table) <- c("Barcode", "UMI", "WT", "MUT", "Body")

sample_names <- c("WT_M_E18_WB1", "WT_M_E18_WB2", "MUT_M_E18_WB1", "MUT_M_E18_WB2", "WT_M_P30_CORT1", "WT_M_P30_CORT2","MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P60_CORT1", "WT_M_P60_CORT2", "MUT_M_P60_CORT1", "MUT_M_P60_CORT2", "WT_M_P120_CORT1", "WT_M_P120_CORT2", "MUT_M_P120_CORT1","MUT_M_P120_CORT2", "MUT_F_P30_CORT1", "MUT_F_P30_CORT2", "MUT_F_P60_CORT1", "MUT_F_P60_CORT2", "MUT_F_P150_CORT1", "MUT_F_P150_CORT2", "MUT_F_P150_CORT3", "MUT_F_P150_CORT4", "WT_F_P30_CORT1", "WT_F_P30_CORT2", "WT_F_P60_CORT1", "WT_F_P60_CORT2", "WT_F_P150_CORT1", "WT_F_P150_CORT2", "WT_F_P150_CORT3", "WT_F_P150_CORT4", "MUT_F_E18_WB1", "MUT_F_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2")
#parse by sex
all.cortex.combined$Sex <- plyr::mapvalues(
  x = all.cortex.combined$orig.ident, 
  from = c(sample_names), 
  to = c("Male", "Male", "Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female")
)
#parse by condition
all.cortex.combined$Condition <- plyr::mapvalues(
  x = all.cortex.combined$orig.ident, 
  from = c(sample_names), 
  to = c("WT", "WT", "MUTANT", "MUTANT", "WT", "WT", "MUTANT", "MUTANT","WT", "WT", "MUTANT", "MUTANT","WT", "WT", "MUTANT", "MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","WT","WT","WT","WT","WT","WT","WT","WT","MUTANT", "MUTANT", "WT","WT")
)

#parse by Age
all.cortex.combined$Age <- plyr::mapvalues(
  x = all.cortex.combined$orig.ident, 
  from = c(sample_names), 
  to = c("E18", "E18","E18","E18","P30","P30","P30","P30","P60","P60","P60","P60","P120","P120","P120","P120","P30", "P30", "P60", "P60", "P150", "P150","P150", "P150","P30", "P30", "P60", "P60", "P150", "P150","P150", "P150","E18", "E18","E18","E18")
)

Idents(all.cortex.combined) <- "celltype.call"
levels(all.cortex.combined) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")

## Adding in Mecp2 transcript expression information into Seurat Object in counts data and in meta data ##

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "dwtools", "devtools")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

#headers for the txt file
names(counts_table) <- c("Barcode", "UMI", "WT", "MUT", "Body")

Mecp2_wt_mut_counts = counts_table
# Adding body counts into WT or MUT counts
new.counts.wt = sapply(1:nrow(Mecp2_wt_mut_counts), function(x) {
  ifelse(Mecp2_wt_mut_counts$Body[x]>0 & Mecp2_wt_mut_counts$WT[x]>0, 
         Mecp2_wt_mut_counts$Body[x]+Mecp2_wt_mut_counts$WT[x],
         Mecp2_wt_mut_counts$WT[x])
})

new.counts.mut = sapply(1:nrow(Mecp2_wt_mut_counts), function(x) {
  ifelse(Mecp2_wt_mut_counts$Body[x]>0 & Mecp2_wt_mut_counts$MUT[x]>0, 
         Mecp2_wt_mut_counts$Body[x]+Mecp2_wt_mut_counts$MUT[x],
         Mecp2_wt_mut_counts$MUT[x])
})

Mecp2_wt_mut_counts$WT = new.counts.wt
Mecp2_wt_mut_counts$MUT = new.counts.mut

# checking to make sure no duplicate barcodes - answer should be true
is.unique(Mecp2_wt_mut_counts$Barcode)

# making a data frame with just one vector representing all barcodes present in Seurat object by removing everything after the - symbol in the barcode names in the Seurat object
dataFrame = data.frame(Barcode = sub("-.*", "", names(Idents(all.cortex.combined))))

# adding in an id column to make sure we can sort the data back to the original order found in the Seurat Object
dataFrame$id = 1:nrow(dataFrame)

#checking to make sure all barcodes in Seurat Object are unique
is.unique(dataFrame$Barcode)

# merging Mecp2 counts with Barcodes in the order that they appear in the Seurat Object
merged = merge(dataFrame, Mecp2_wt_mut_counts, by="Barcode", all=T)

# reordering to original Seurat Object order
merged = merged[order(merged$id),]

# converting NAs to zeros
merged[is.na(merged)] <- 0

# Number of barcodes that are not in Seurat object
number = nrow(merged) - all.cortex.combined@assays$RNA@counts@Dim[2] # x number of Barcodes in the Mecp2 WT vs. MUT counts file did not have corresponding Barcodes in the Seurat Object

# Removing rows with barcodes that are not in Seurat object
merged = merged[c(1:(nrow(merged)-number)),]

nrow(merged) # should equal length of all metadata 

### Adding counts to Seurat metadata 

all.cortex.combined$WT_Mecp2 = merged$WT

all.cortex.combined$MUT_Mecp2 = merged$MUT

### Adding counts to Seurat counts as two new genes

all.cortex.combined@assays$RNA@counts = rbind(all.cortex.combined@assays$RNA@counts, merged$WT, merged$MUT)
nrow(all.cortex.combined@assays$RNA@counts)
rownames(all.cortex.combined@assays$RNA@counts) = c(rownames(all.cortex.combined@assays$RNA@counts)[-c(19665:19666)], "Mecp2-WT", "Mecp2-MUT")

# checking to make sure they are added
all.cortex.combined@assays$RNA@counts[19665:19666, 1:5]

#E18 <- subset(x = all.cortex.combined, subset = orig.ident == c("MUT_F_E18_WB1", "MUT_F_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2"))
all.cortex.combined <- NormalizeData(all.cortex.combined)
all.cortex.combined <- ScaleData(all.cortex.combined)
all.cortex.combined <- FindVariableFeatures(all.cortex.combined, selection.method = "vst", nfeatures = 3000)
all.cortex.combined <- RunPCA(object = all.cortex.combined, verbose = FALSE)
all.cortex.combined <- RunUMAP(object = all.cortex.combined, dims = 1:20, verbose = FALSE)
all.cortex.combined <- FindNeighbors(object = all.cortex.combined, dims = 1:20, verbose = FALSE)
all.cortex.combined <- FindClusters(object = all.cortex.combined, verbose = FALSE)
DimPlot(object = all.cortex.combined, label = TRUE, group.by = 'celltype.call') + NoLegend() + ggtitle("sctransform")# saving new Seurat object
FeaturePlot_scCustom(seurat_object = all.cortex.combined, features = 'WT_Mecp2')
FeaturePlot_scCustom(seurat_object = all.cortex.combined, features = 'MUT_Mecp2')
FeaturePlot_scCustom(seurat_object = all.cortex.combined, features = 'Mecp2', split.by = "Sex")

#Function counts the percent of total cells that express specific genes
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

PrctCellExpringGene(all.cortex.combined, genes =c("WT-Mecp2","MUT-Mecp2"), group.by = "all")

GOI1 <- 'WT_Mecp2' #you will have to name your first gene here, im choosing PDX1 as an example
GOI2 <- 'MUT_Mecp2' #you will have to name your second gene here, im choosing INS as an example
GOI1.cutoff <- 1 #Assumption: gene count cutoff is 1, assuming atleast 1 count is REAL
GOI2.cutoff <- 1 #Assumption: gene count cutoff is 1, assuming atleast 1 count is REAL

# Time to party
GOI1.cells <- length(which(FetchData(E18, vars = GOI1) > GOI1.cutoff))
GOI2.cells <- length(which(FetchData(E18, vars = GOI2) > GOI2.cutoff))
GOI1_GOI2.cells <- length(which(FetchData(E18, vars = GOI2) > GOI2.cutoff & FetchData(E18, vars = GOI1) > GOI1.cutoff))
all.cells.incluster <- table(E18$orig.ident)
GOI1.cells/all.cells.incluster*100 # Percentage of cells in Beta that express GOI1
GOI2.cells/all.cells.incluster*100 #Percentage of cells in Beta that express GOI2
GOI1_GOI2.cells/all.cells.incluster*100 #Percentage of cells in Beta that co-express GOI1 + GOI2

#Marker Genes
cell_markers_manual <- c("Plch2","Sst","Vip", "Pvalb", "Slc17a8", "Macc1", "Rorb", "Fezf2", "Rprm", "Aqp4", "Rassf10", "Kcnj8", "Slc17a7", "Gad2", "Aspa")

FeaturePlot_scCustom(seurat_object = all.cortex.combined, features = cell_markers_manual)
DotPlot_scCustom(seurat_object = all.cortex.combined, features = cell_markers_manual, group.by = "celltype.call", x_lab_rotate = TRUE, flip_axes = TRUE)
ggsave("/Users/osman/Desktop/LaSalle_lab/Seurat_figures/celltype_markers.pdf")
markers_df <- FindAllMarkers(
  object = all.cortex.combined, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  thresh.use = 0.25
)
markers_all_single <- markers_df[markers_df$gene %in% names(table(markers_df$gene))[table(markers_df$gene) == 1],]
top5 <- markers_all_single %>% group_by(cluster) %>% top_n(5, avg_log2FC)
dim(top5)
DoHeatmap(
  object = all.cortex.combined, 
  features = top5$gene,
  labels = FALSE) 

save(all.cortex.combined, file="/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all.cortex.combined.RData")
save(all.cortex.combined, file=glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Parsing_Mecp2_trnx_expression/{all.cortex.combined.name}/{all.cortex.combined.name}_with_Mecp2_WT_MUT.RData"))

