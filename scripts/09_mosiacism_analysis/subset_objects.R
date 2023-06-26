library(Seurat)
library(dplyr)
library(scCustomize)
library(patchwork)
library(glue)
library(stringr)

##########################
## load complete object ##
##########################
load('/Users/osman/Desktop/LaSalle_lab/Seurat_objects/all.cortex.combined.RData')
Idents(all.cortex.combined) <- 'celltype.call'
DefaultAssay(all.cortex.combined) <- 'RNA'
levels(all.cortex.combined) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")

########################
## subset adult males ##
########################
all.male.cortex <- subset(x = all.cortex.combined, subset = Sex == 'Male')
all.male.cortex <- subset(x = all.male.cortex, subset = Age != 'E18')

############################
## load adult male counts ##
############################
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
combined_male_counts <- read.csv(glue('{base_path}/combined_male_unique.csv'))
combined_male_counts[1:2] <- NULL
combined_male_counts <- combined_male_counts[!(combined_male_counts$Condition == "Mecp2_MUT" & combined_male_counts$WT_Mecp2 > 0), ]
combined_male_counts <- combined_male_counts[!(combined_male_counts$Condition == "Mecp2_WT" & combined_male_counts$MUT_Mecp2 > 0), ]
# Create the 'metadata' column
combined_male_counts$metadata <- paste(combined_male_counts$Timepoint,
                                       combined_male_counts$Barcode,
                                       combined_male_counts$Condition,
                                       combined_male_counts$Sex,
                                       combined_male_counts$Timepoint,
                                       combined_male_counts$Sample,
                                       sep = "_")
# Update the 'metadata' column
combined_male_counts$metadata <- gsub("_Mecp2_", "-", combined_male_counts$metadata)

##########################################
## add counts to seurat object metadata ##
##########################################
# Extract the barcodes from Seurat object's meta.data
#barcodes <- str_extract(rownames(all.male.cortex@meta.data), "(?<=_).*(?=-)")
barcodes <- rownames(all.male.cortex@meta.data)
# Create a new column 'Mecp2_allele' in Seurat object's meta.data
all.male.cortex@meta.data$Mecp2_allele <- ""
# Find matching barcodes and transfer corresponding info
matching_rows <- combined_male_counts$metadata %in% barcodes
# Find the non-matching items
non_matching_rows <- combined_male_counts$metadata[!(combined_male_counts$metadata %in% barcodes)]
all.male.cortex@meta.data$Mecp2_allele[matching_rows] <- combined_male_counts$Condition[match(barcodes[matching_rows], combined_male_counts$metadata)]
all.male.cortex@meta.data$Mecp2_allele <- ifelse(is.na(na_if(all.male.cortex@meta.data$Mecp2_allele, "N/A")), "", all.male.cortex@meta.data$Mecp2_allele)
Idents(all.male.cortex) <- 'Mecp2_allele'
all.male.cortex <- subset(x = all.male.cortex, idents = c("Mecp2_WT", "Mecp2_MUT"))
Idents(all.male.cortex) <- 'celltype.call'
levels(all.male.cortex) <- c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5","Peri", "Endo", "Oligo","Astro","Non-neuronal")
DimPlot_scCustom(seurat_object = all.male.cortex, pt.size = 0.6, split.by = 'Condition', group.by = 'Mecp2_allele', order = TRUE, label.box = TRUE, colors_use = c('green3', 'red'))
ggplot2::ggsave(glue("{base_path}/allmale_Mecp2allele_green_UMAP.tiff"),
                device = NULL,
                height = 10,
                width = 12)
write.csv(combined_male_counts, file = glue('{base_path}/adult_male_counts.csv'), row.names = TRUE)
save(all.male.cortex, file = glue('{base_path}/all.male.cortex.parsed.RData'))