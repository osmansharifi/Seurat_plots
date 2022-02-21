## Adding in Mecp2 transcript expression information into Seurat Object in counts data and in meta data ##

# install dwtools package
library(devtools)
install_github("jangorecki/dwtools")

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat", "dwtools")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

s.obj.name = "all_female_P150" #change this to the name of the Seurat object you're working with

load(glue::glue("/Users/karineier/Documents/Mecp2/scRNA-seq/{s.obj.name}.RData"))

s.obj = all_female_P150 #change this to the name of the Seurat object you're working with

# reading in count data for Mecp2 transcript expression wt vs.  mut
Mecp2_wt_mut_counts = read.table(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Parsing_Mecp2_trnx_expression/{s.obj.name}/counts_table.txt"), header=T)

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

s.obj$WT_Mecp2_counts = merged$WT

s.obj$MUT_Mecp2_counts = merged$MUT

### Adding counts to Seurat counts as two new genes

s.obj@assays$RNA@counts = rbind(s.obj@assays$RNA@counts, merged$WT, merged$MUT)
nrow(s.obj@assays$RNA@counts)
rownames(s.obj@assays$RNA@counts) = c(rownames(s.obj@assays$RNA@counts)[-c(7531:7532)], "Mecp2-WT", "Mecp2-MUT")

# checking to make sure they are added
s.obj@assays$RNA@counts[7530:7532, 1:5]

# saving new Seurat object
save(s.obj, file=glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Parsing_Mecp2_trnx_expression/{s.obj.name}/{s.obj.name}_with_Mecp2_WT_MUT.RData"))

