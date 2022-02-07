### Differential Expression Analysis of scRNA-seq Data - 00 - Splitting Seurat Object ###

# Splitting Seurat Object by Activated and Unactivated Neurons

# set up 

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

s.obj.name = "rett_P60_with_labels_proportions" #change this to the name of the Seurat object you're working with

load(glue::glue("/Users/karineier/Documents/Mecp2/scRNA-seq/{s.obj.name}.rda"))

s.obj = experiment.aggregate #change this to the name of the Seurat object you're working with

dir.create(glue::glue("{s.obj.name}"))

setwd(glue::glue("{s.obj.name}"))

### Split Seurat object into activated and unactivated neurons ###

activated.neurons.df = ifelse((s.obj$celltype.call == "L2_3_IT" | 
                                 s.obj$celltype.call == "L4" | 
                                 s.obj$celltype.call == "L5" | 
                                 s.obj$celltype.call == "L6" |
                                 s.obj$celltype.call ==  "Pvalb" | 
                                 s.obj$celltype.call == "Lamp5" | 
                                 s.obj$celltype.call == "Vip" | 
                                 s.obj$celltype.call == "Sst" | 
                                 s.obj$celltype.call == "Sncg") &
                                s.obj@assays$RNA@counts["Fos",]>0, "-activated", "")

activated.names = sapply(1:length(s.obj$celltype.call), function(x){
  new.id = gsub("$", glue::glue("{activated.neurons.df[[x]]}"), s.obj$celltype.call[[x]])
  return(new.id)
})

s.obj$cell.type = activated.names

Idents(s.obj) = "cell.type"

save(s.obj, activated.neurons.df, activated.names, file="DEanalysis_00.RData")
