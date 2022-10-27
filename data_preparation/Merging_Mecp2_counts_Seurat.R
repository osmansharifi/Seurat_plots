## Adding in Mecp2 transcript expression information into Seurat Object in counts data and in meta data ##

# install dwtools package
#install_github("jangorecki/dwtools")

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "dwtools", "devtools")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

s.obj.name = "all_embryonic" #change this to the name of the Seurat object you're working with

load(glue::glue("/Users/osman/Desktop/All_female_samples/all_female_cortex_labeled.RData"))

s.obj = mousecortexsca #change this to the name of the Seurat object you're working with
s.obj = all_female.query
# reading in count data for Mecp2 transcript expression wt vs.  mut
counts_table <- read.table("/Users/osman/Desktop/LaSalle_lab/Scripts/E18_script/E18_Male_cortex/master_E18_Mecp2_counts.txt", sep="\t", header=FALSE)
names(counts_table) <- c("Barcode", "UMI", "WT", "MUT", "Body")
#Mecp2_wt_mut_counts = read.table(glue::glue("/Users/osman/Downloads/E18_mut.alleler"), header=T)
Mecp2_wt_mut_counts = counts_table
Mecp2_wt_mut_counts = Mecp2
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
dataFrame = data.frame(Barcode = sub("-.*", "", names(Idents(s.obj))))

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
number = nrow(merged) - s.obj@assays$RNA@counts@Dim[2] # 12 Barcodes in the Mecp2 WT vs. MUT counts file did not have corresponding Barcodes in the Seurat Object

# Removing rows with barcodes that are not in Seurat object
merged = merged[c(1:(nrow(merged)-number)),]

nrow(merged) # should equal length of all metadata 

### Adding counts to Seurat metadata 

s.obj$WT_Mecp2 = merged$WT

s.obj$MUT_Mecp2 = merged$MUT

### Adding counts to Seurat counts as two new genes

s.obj@assays$RNA@counts = rbind(s.obj@assays$RNA@counts, merged$WT, merged$MUT)
nrow(s.obj@assays$RNA@counts)
rownames(s.obj@assays$RNA@counts) = c(rownames(s.obj@assays$RNA@counts)[-c(19667:19668)], "Mecp2-WT", "Mecp2-MUT")

# checking to make sure they are added
s.obj@assays$RNA@counts[7530:7532, 1:5]

E18 <- subset(x = s.obj, subset = orig.ident == c("MUT_F_E18_WB1", "MUT_F_E18_WB2", "WT_F_E18_WB1", "WT_F_E18_WB2"))
E18 <- RunPCA(object = E18, verbose = FALSE)
E18 <- RunUMAP(object = E18, dims = 1:20, verbose = FALSE)
E18 <- FindNeighbors(object = E18, dims = 1:20, verbose = FALSE)
E18 <- FindClusters(object = E18, verbose = FALSE)
DimPlot(object = E18, label = TRUE, group.by = 'predicted.id') + NoLegend() + ggtitle("sctransform")# saving new Seurat object
FeaturePlot_scCustom(seurat_object = E18, features = 'WT_Mecp2')
FeaturePlot_scCustom(seurat_object = E18, features = 'MUT_Mecp2')
FeaturePlot_scCustom(seurat_object = E18, features = 'Mecp2')

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

save(s.obj, file=glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Parsing_Mecp2_trnx_expression/{s.obj.name}/{s.obj.name}_with_Mecp2_WT_MUT.RData"))


