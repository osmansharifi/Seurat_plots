library(Seurat)

#Load data
load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all_rett_mouse_cortex.RData")
counts_table <- read.table("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/data_preparation/master_sequence_alleler.txt", sep="\t", header=FALSE)

sample_names <- c("WT_M_E18_WB1", "WT_M_E18_WB2", "MUT_M_E18_WB1", "MUT_M_E18_WB2", "WT_M_P30_CORT1", "WT_M_P30_CORT2","MUT_M_P30_CORT1", "MUT_M_P30_CORT2", "WT_M_P60_CORT1", "WT_M_P60_CORT2", "MUT_M_P60_CORT1", "MUT_M_P60_CORT2", "WT_M_P120_CORT1", "WT_M_P120_CORT2", "MUT_M_P120_CORT1","MUT_M_P120_CORT2", "MUT_F_P30_CORT1", "MUT_F_P30_CORT2", "MUT_F_P60_CORT1", "MUT_F_P60_CORT2", "MUT_F_P150_CORT1", "MUT_F_P150_CORT2", "MUT_F_P150_CORT3", "MUT_F_P150_CORT4", "WT_F_P30_CORT1", "WT_F_P30_CORT2", "WT_F_P60_CORT1", "WT_F_P60_CORT2", "WT_F_P150_CORT1", "WT_F_P150_CORT2", "WT_F_P150_CORT3", "WT_F_P150_CORT4", "MUT_F_E18_WB1", "MUT_F_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2")
#parse by sex
all_rett_mouse_cortex$Sex <- plyr::mapvalues(
  x = all_rett_mouse_cortex$orig.ident, 
  from = c(sample_names), 
  to = c("Male", "Male", "Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Male","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female")
)
#parse by condition
all_rett_mouse_cortex$Condition <- plyr::mapvalues(
  x = all_rett_mouse_cortex$orig.ident, 
  from = c(sample_names), 
  to = c("WT", "WT", "MUTANT", "MUTANT", "WT", "WT", "MUTANT", "MUTANT","WT", "WT", "MUTANT", "MUTANT","WT", "WT", "MUTANT", "MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","WT","WT","WT","WT","WT","WT","WT","WT","MUTANT", "MUTANT", "WT","WT")
)

#parse by Age
all_rett_mouse_cortex$Age <- plyr::mapvalues(
  x = all_rett_mouse_cortex$orig.ident, 
  from = c(sample_names), 
  to = c("E18", "E18","E18","E18","P30","P30","P30","P30","P60","P60","P60","P60","P120","P120","P120","P120","P30", "P30", "P60", "P60", "P150", "P150","P150", "P150","P30", "P30", "P60", "P60", "P150", "P150","P150", "P150","E18", "E18","E18","E18")
)

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
dataFrame = data.frame(Barcode = sub("-.*", "", names(Idents(all_rett_mouse_cortex))))

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
number = nrow(merged) - all_rett_mouse_cortex@assays$RNA@counts@Dim[2] # x number of Barcodes in the Mecp2 WT vs. MUT counts file did not have corresponding Barcodes in the Seurat Object

# Removing rows with barcodes that are not in Seurat object
merged = merged[c(1:(nrow(merged)-number)),]

nrow(merged) # should equal length of all metadata 

### Adding counts to Seurat metadata 

all_rett_mouse_cortex$WT_Mecp2 = merged$WT

all_rett_mouse_cortex$MUT_Mecp2 = merged$MUT

### Adding counts to Seurat counts as two new genes

all_rett_mouse_cortex@assays$RNA@counts = rbind(all_rett_mouse_cortex@assays$RNA@counts, merged$WT, merged$MUT)
nrow(all_rett_mouse_cortex@assays$RNA@counts)
rownames(all_rett_mouse_cortex@assays$RNA@counts) = c(rownames(all_rett_mouse_cortex@assays$RNA@counts)[-c(19667:19668)], "Mecp2-WT", "Mecp2-MUT")

# checking to make sure they are added
all_rett_mouse_cortex@assays$RNA@counts[19667:19668, 1:5]

#E18 <- subset(x = all_rett_mouse_cortex, subset = orig.ident == c("MUT_F_E18_WB1", "MUT_F_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2"))
all_rett_mouse_cortex <- NormalizeData(all_rett_mouse_cortex)
all_rett_mouse_cortex <- ScaleData(all_rett_mouse_cortex)
all_rett_mouse_cortex <- FindVariableFeatures(all_rett_mouse_cortex, selection.method = "vst", nfeatures = 3000)
all_rett_mouse_cortex <- RunPCA(object = all_rett_mouse_cortex, verbose = FALSE)
all_rett_mouse_cortex <- RunUMAP(object = all_rett_mouse_cortex, dims = 1:20, verbose = FALSE)
all_rett_mouse_cortex <- FindNeighbors(object = all_rett_mouse_cortex, dims = 1:20, verbose = FALSE)
all_rett_mouse_cortex <- FindClusters(object = all_rett_mouse_cortex, verbose = FALSE)
DimPlot(object = all_rett_mouse_cortex, label = TRUE, group.by = 'celltype.call') + NoLegend() + ggtitle("sctransform")# saving new Seurat object
FeaturePlot_scCustom(seurat_object = all_rett_mouse_cortex, features = 'WT_Mecp2')
FeaturePlot_scCustom(seurat_object = all_rett_mouse_cortex, features = 'MUT_Mecp2')
FeaturePlot_scCustom(seurat_object = all_rett_mouse_cortex, features = 'Mecp2')

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

PrctCellExpringGene(E18, genes =c("Mecp2","AC149090.1"), group.by = "all")

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

save(all_rett_mouse_cortex, file=glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Parsing_Mecp2_trnx_expression/{all_rett_mouse_cortex.name}/{all_rett_mouse_cortex.name}_with_Mecp2_WT_MUT.RData"))




save(all_rett_mouse_cortex, file="/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all_rett_mouse_cortex.RData")
