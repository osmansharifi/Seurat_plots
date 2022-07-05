library(Seurat)
library(scCustomize)
library(plyr)

load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all_male_mouse_cortex.RData")
load("/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all_female_mouse_cortex.RData")
mouse_rettcort <- merge(all_male, y = all_female.query, add.cell.ids = c("male", "female"), project = "mouse_rettcort_all",
                         merge.data = TRUE)

samples <- c("MUT_F_E18_WB1","MUT_F_E18_WB2", "MUT_F_P150_CORT1","MUT_F_P150_CORT2","MUT_F_P150_CORT3","MUT_F_P150_CORT4","MUT_F_P30_CORT1","MUT_F_P30_CORT2","MUT_F_P60_CORT1","MUT_F_P60_CORT2","MUT_M_E18_WB1","MUT_M_E18_WB2","MUT_M_P120_CORT1","MUT_M_P120_CORT2","MUT_M_P30_CORT1","MUT_M_P30_CORT2","MUT_M_P60_CORT1","MUT_M_P60_CORT2","WT_F_E18_WB1","WT_F_E18_WB2","WT_F_P150_CORT1","WT_F_P150_CORT2","WT_F_P150_CORT3","WT_F_P150_CORT4","WT_F_P30_CORT1","WT_F_P30_CORT2","WT_F_P60_CORT1","WT_F_P60_CORT2","WT_M_E18_WB1","WT_M_E18_WB2","WT_M_P120_CORT1","WT_M_P120_CORT2","WT_M_P30_CORT1","WT_M_P30_CORT2","WT_M_P60_CORT1","WT_M_P60_CORT2" )

sex <- c("Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Male", "Male", "Male","Male","Male","Male","Male", "Male","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Male","Male","Male","Male","Male","Male","Male","Male")

condition <- c("MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","MUTANT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT")

age <- c("E18", "E18", "P150", "P150", "P150", "P150", "P30", "P30", "P60", "P60", "E18", "E18", "P120", "P120", "P30", "P30", "P60", "P60", "E18", "E18", "P150", "P150", "P150", "P150","P30", "P30", "P60", "P60", "E18", "E18", "P120", "P120", "P30", "P30", "P60", "P60")

sactime <- c("Embryonic", "Embryonic", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Embryonic", "Embryonic", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Embryonic", "Embryonic", "Post-natal", "Post-natal", "Post-natal", "Post-natal","Post-natal", "Post-natal", "Post-natal", "Post-natal", "Embryonic", "Embryonic", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal", "Post-natal")

#parse by sex
mouse_rettcort$Sex <- plyr::mapvalues(
  x = mouse_rettcort$orig.ident, 
  from = c(samples), 
  to = c(sex)
)
#parse by condition
mouse_rettcort$Condition <- plyr::mapvalues(
  x = mouse_rettcort$orig.ident, 
  from = c(samples), 
  to = c(condition)
)

#parse by Age
mouse_rettcort$Age <- plyr::mapvalues(
  x = mouse_rettcort$orig.ident, 
  from = c(samples), 
  to = c(age)
)

#parse by sactime
mouse_rettcort$SacTime <- plyr::mapvalues(
  x = mouse_rettcort$orig.ident, 
  from = c(samples), 
  to = c(sactime)
)

#Split object into two
obj.list <- SplitObject(mouse_rettcort, split.by = "SacTime")
mouse_postnatal_cortex <- obj.list$`Post-natal`
mouse_embryonic_cortex <- obj.list$Embryonic

save(mouse_embryonic_cortex, file="mouse_embryonic_cortex.RData")
save(mouse_postnatal_cortex, file="mouse_postnatal_cortex.RData")

table(mouse_rettcort$orig.ident)
LoadData("humancortexref")
#The RunAzimuth function can take a Seurat object as input
mouse_rettcort_test <- RunAzimuth(mouse_postnatal_cortex, reference = "mousecortexref")
#Normalize data
all_male <- NormalizeData(all_male)
all_male <- ScaleData(all_male)
Idents(all_male) <- "predicted.subclass"

DimPlot_scCustom(seurat_object = all_male)
