### Differential Expression Analysis of scRNA-seq Data - 00 - Splitting Seurat Object ###

# Splitting Seurat Object by Activated and Unactivated Neurons

# set up 

packages <- c("tidyr", "openxlsx", "glue", "magrittr", "Seurat")
stopifnot(suppressMessages(sapply(packages, require, character.only=TRUE)))

s.obj.name = "all.cortex.combined" #change this to the name of the Seurat object you're working with

load(glue::glue("/Users/karineier/Documents/Mecp2/scRNA-seq/{s.obj.name}.RData"))

s.obj = all.cortex.combined #change this to the name of the Seurat object you're working with

dir.create(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/Mecp2e1_parsed"))

setwd(glue::glue("/Users/karineier/Documents/GitHub/snRNA-seq-pipeline/Differential_expression/Mecp2e1_parsed"))

# checking for Mecp2-WT and Mecp2-MUT data 

sum(s.obj@assays$RNA@counts[which(grepl("Mecp2", rownames(s.obj@assays$RNA@counts))),])

sum(s.obj$WT_Mecp2)
sum(s.obj$MUT_Mecp2)

length(which(s.obj@assays$RNA@counts[which(grepl("Mecp2", rownames(s.obj@assays$RNA@counts))),]==0 & s.obj$Mecp2_expr != "NONE"))

print(paste("The number of parsed Mecp2 transcripts was", sum(s.obj$WT_Mecp2)+sum(s.obj$MUT_Mecp2), "which was", ((sum(s.obj$WT_Mecp2)+sum(s.obj$MUT_Mecp2))/sum(s.obj@assays$RNA@counts[which(grepl("Mecp2", rownames(s.obj@assays$RNA@counts))),]))*100, "% of all Mecp2 transcripts in the Seurat object", sep=" "))

### Split Seurat object into Mecp2 expressing cells ###

mecp2.expression.df = ifelse(s.obj$MUT_Mecp2>0, "MUTANT", ifelse(s.obj$WT_Mecp2>0 & s.obj$MUT_Mecp2==0, "WT", ifelse(s.obj$MUT_Mecp2==0 & s.obj$WT_Mecp2==0, "NONE", "BOTH")))
mecp2.expression.mut = ifelse(s.obj$MUT_Mecp2>0, "MUTANT", "NOT")
mecp2.expression.wt = ifelse(s.obj$WT_Mecp2>0, "WT", "NOT")

length(which(mecp2.expression.df=="MUTANT"))
length(which(mecp2.expression.mut=="MUTANT"))
length(which(mecp2.expression.df=="WT"))
length(which(mecp2.expression.wt=="WT"))

length(which(mecp2.expression.df=="MUTANT" & mecp2.expression.mut=="MUTANT"))
length(which(mecp2.expression.df!="WT" & mecp2.expression.wt=="WT"))

cells.mut <- length(mecp2.expression.df[which(mecp2.expression.df=="MUTANT")])

print(glue::glue("The number of cells expressing mutant Mecp2-e1 is {cells.mut}"))

cells.wt <- length(mecp2.expression.df[which(mecp2.expression.df=="WT")])

print(glue::glue("The number of cells expressing wild-type Mecp2-e1 is {cells.wt}"))

cells.both <- length(mecp2.expression.df[which(mecp2.expression.df=="BOTH")])

print(glue::glue("The number of cells expressing both mutant and wild-type Mecp2-e1 is {cells.both}"))
                                 
cells.none <- length(mecp2.expression.df[which(mecp2.expression.df=="NONE")])

print(glue::glue("The number of cells expressing no Mecp2-e1 is {cells.none}"))

# checking to make sure everything adds up - these next two lines should have the same result

ncol(s.obj@assays$RNA@counts)

cells.mut + cells.wt + cells.none

# adding new column in metadata describing Mecp2 expression

s.obj$Mecp2_expr <- mecp2.expression.df

## creating expression and design matrices

# filter out cells with no Mecp2 expression and filter out cells from E18 mice
expr_matrix = s.obj@assays$RNA@counts[,which(s.obj$Mecp2_expr != "NONE")]


design = data.frame(cell_type = s.obj$celltype.call,
                    sample_ID = s.obj$orig.ident,
                    cell_cycle = s.obj$cell.cycle,
                    cell_ID = colnames(s.obj@assays$RNA@counts),
                    genotype = s.obj$Condition,
                    percent.mito = s.obj$percent.mito,
                    age = s.obj$Age,
                    Mecp2_expressing = ifelse(s.obj$Mecp2_expr != "NONE", "YES", "NO"),
                    Mecp2e1_expression = s.obj$Mecp2_expr, 
                    sex = s.obj$Sex,
                    Mecp2e1_MUT = s.obj$MUT_Mecp2,
                    Mecp2e1_WT = s.obj$WT_Mecp2
                    )

design = design %>%
  dplyr::mutate_if(is.character, as.factor)

design$genotype = factor(design$genotype, levels=c("WT", "MUTANT"))

design$age = factor(design$age, levels=c("E18", "P30", "P60", "P120", "P150"))

design$Mecp2e1_expression = factor(design$Mecp2e1_expression, levels=c("WT", "MUTANT", "NONE"))

design = design %>%
  dplyr::filter(Mecp2_expressing == "YES") 

design$Mecp2e1_expression = droplevels(design$Mecp2e1_expression)

save(s.obj, file = "DEanalysis_Mecp2_parsing_Seurat_Object_00.RData")

save(expr_matrix, design, file="DEanalysis_Mecp2_parsing_00.RData")
