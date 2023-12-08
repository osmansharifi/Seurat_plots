library(dplyr)
library(GeneOverlap)

#####################
## Load mouse DEGs ##
#####################
mouse_DEGs <- read.csv('/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/total_sig_mouse_DEGs_limmaVoom.csv')
mouse_DEGs$X <- NULL
head(mouse_DEGs)
dim(mouse_DEGs)
mouse_male_DEGs <- mouse_DEGs[mouse_DEGs$Sex !='females', ]
dim(mouse_male_DEGs)
# Filter rows where 'adj.P.Val' is 0.05 or lower
mouse_male_DEGs <- mouse_male_DEGs[mouse_male_DEGs$adj.P.Val <= 0.05, ]
# Convert 'SYMBOL' column to uppercase
mouse_male_DEGs$SYMBOL <- toupper(mouse_male_DEGs$SYMBOL)
# Create the 3 broad cell classes
mouse_male_DEGs <- mouse_male_DEGs %>%
  mutate(broad_class = case_when(
    Cell_Type %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Vip") ~ "GABAergic",
    Cell_Type %in% c("L2_3_IT", "L4", "L5", "L6") ~ "Glutamatergic",
    Cell_Type %in% c("Astro", "Non-neuronal", "Oligo") ~ "Non-neuronal",
    TRUE ~ "Other"
  ))

#####################
## Load human DEGs ##
#####################
human_DEGs <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/total_sig_human_DEGs_limmaVoom.csv")
# Filter rows where 'adj.P.Val' is 0.05 or lower
human_female_DEGs <- human_DEGs[human_DEGs$adj.P.Val <= 0.05, ]
# Create the 3 broad cell classes
human_female_DEGs <- human_female_DEGs %>%
  mutate(broad_class = case_when(
    Cell_Type %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Vip") ~ "GABAergic",
    Cell_Type %in% c("L2_3_IT", "L5_6", "L6") ~ "Glutamatergic",
    Cell_Type %in% c("Astro", "OPC", "Oligo", "Non-neuronal") ~ "Non-neuronal",
    TRUE ~ "Other"
  ))
head(human_female_DEGs)

###################################################
## Overlap the human DEGs with female mouse DEGs ##
###################################################
# Extract significant DEGs
Glut_human_female_DEGs <- human_female_DEGs$SYMBOL[human_female_DEGs$broad_class == 'Glutamatergic']
Glut_mouse_male_DEGs <- mouse_male_DEGs$SYMBOL[mouse_male_DEGs$broad_class == 'Glutamatergic']
GABA_human_female_DEGs <- human_female_DEGs$SYMBOL[human_female_DEGs$broad_class == 'GABAergic']
GABA_mouse_male_DEGs <- mouse_male_DEGs$SYMBOL[mouse_male_DEGs$broad_class == 'GABAergic']
Non_neuronal_human_female_DEGs <- human_female_DEGs$SYMBOL[human_female_DEGs$broad_class == 'Non-neuronal']
Non_neuronal_mouse_male_DEGs <- mouse_male_DEGs$SYMBOL[mouse_male_DEGs$broad_class == 'Non-neuronal']
intersection_Glut <- intersect(Glut_human_female_DEGs,Glut_mouse_male_DEGs)
intersection_GABA <- intersect(GABA_human_female_DEGs,GABA_mouse_male_DEGs)
intersection_non_neuronal <- intersect(Non_neuronal_human_female_DEGs,Non_neuronal_mouse_male_DEGs)
# Create a Venn diagram
pdf(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/male_venn_Glutamatergic.pdf"))
temp <- venn.diagram(
  x = list(
    Glut_human = Glut_human_female_DEGs,
    Glut_mouse = Glut_mouse_male_DEGs
  ),
  category.names = c("Female human Glutamatergic DEGs", "Male mouse Glutamatergic DEGs"),
  main = 'Glutamatergic DEGs from human and mouse ',
  #filename = glue("{base_path}/broad_group_analysis/venn_glutamatergic.pdf"),
  filename = NULL,
  col = c('#F79120', '#372A82'), 
  fill = c('#F79119', '#372A81'),
  cat.cex = 1.2,
  cat.fontface = "bold",
  euler.d = TRUE,
  disable.logging = TRUE,
  hyper.test = TRUE
)
grid.draw(temp)
dev.off()

#GABAergic cells
pdf(glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/male_venn_GABAergic.pdf"))
temp <- venn.diagram(
  x = list(
    GABA_human = GABA_human_female_DEGs,
    GABA_mouse = GABA_mouse_male_DEGs
  ),
  category.names = c("Female human GABAergic DEGs", "Male mouse GABAergic DEGs"),
  main = 'GABAergic DEGs from human and mouse ',
  #filename = glue("{base_path}/broad_group_analysis/venn_glutamatergic.pdf"),
  filename = NULL,
  col = c('#F79120', '#372A82'), 
  fill = c('#F79119', '#372A81'),
  cat.cex = 1.2,
  cat.fontface = "bold",
  euler.d = TRUE,
  disable.logging = TRUE,
  hyper.test = TRUE
)
grid.draw(temp)
dev.off()

# Make all vectors the same length (pad with NA)
max_length <- max(length(Glut_human_female_DEGs), length(Glut_mouse_male_DEGs))
Glut_human_female_DEGs_overlap <- c(Glut_human_female_DEGs, rep(NA, max_length - length(Glut_human_female_DEGs)))
Glut_mouse_male_DEGs_overlap <- c(Glut_mouse_male_DEGs, rep(NA, max_length - length(Glut_mouse_male_DEGs)))


# Combine into a dataframe
combined_df <- data.frame(Glut_human_female_DEGs_overlap, Glut_mouse_male_DEGs_overlap)

go.obj <- newGeneOverlap(combined_df$Glut_human_female_DEGs_overlap,
                         combined_df$Glut_mouse_male_DEGs_overlap,
                         genome.size = 6500)
go.obj <- testGeneOverlap(go.obj)
getPval(go.obj)
getOddsRatio(go.obj)
getJaccard(go.obj)
getContbl(go.obj)
print(go.obj)
write.csv(as.data.frame(go.obj), file = ('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/geneoverlap_glut.txt'))

base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
enrichR:::.onAttach()
source(glue('{base_path}/GO_ploting_functions.R'))
tryCatch({
  intersection_Glut %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Molecular_Function_2023",
                       "GO_Cellular_Component_2023",
                       "KEGG_2021_Human",
                       "KEGG_2019_Mouse",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
    openxlsx::write.xlsx(file = glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/Glut_intersection_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Mm.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/Glut_intersection_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/Glut_intersection_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10)
},
error = function(error_condition) {
  print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
})

print(glue::glue("The pipeline has finished for samples"))


## GABAergic
# Make all vectors the same length (pad with NA)
max_length <- max(length(GABA_human_female_DEGs), length(GABA_mouse_male_DEGs))
GABA_human_female_DEGs_overlap <- c(GABA_human_female_DEGs, rep(NA, max_length - length(GABA_human_female_DEGs)))
GABA_mouse_male_DEGs_overlap <- c(GABA_mouse_male_DEGs, rep(NA, max_length - length(GABA_mouse_male_DEGs)))


# Combine into a dataframe
combined_df <- data.frame(GABA_human_female_DEGs_overlap, GABA_mouse_male_DEGs_overlap)

go.obj <- newGeneOverlap(combined_df$GABA_human_female_DEGs_overlap,
                         combined_df$GABA_mouse_male_DEGs_overlap,
                         genome.size = 6500)
go.obj <- testGeneOverlap(go.obj)
getPval(go.obj)
getOddsRatio(go.obj)
getJaccard(go.obj)
getContbl(go.obj)
print(go.obj)
write.csv(as.data.frame(go.obj), file = ('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/geneoverlap_glut.txt'))

base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
enrichR:::.onAttach()
source(glue('{base_path}/GO_ploting_functions.R'))
tryCatch({
  intersection_Glut %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Molecular_Function_2023",
                       "GO_Cellular_Component_2023",
                       "KEGG_2021_Human",
                       "KEGG_2019_Mouse",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>%
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis="")) %T>%
    openxlsx::write.xlsx(file = glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/Glut_intersection_enrichr.xlsx")) %>%
    slimGO(tool = "enrichR",
           annoDb = "org.Mm.eg.db",
           plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/Glut_intersection_rrvgo_enrichr.xlsx")) %>%
    GOplot() %>%
    ggplot2::ggsave(glue::glue("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/overlapping_femaleMouse_human/Glut_intersection_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10)
},
error = function(error_condition) {
  print(glue::glue("ERROR: Gene Ontology pipe did not finish for samples"))
})

print(glue::glue("The pipeline has finished for samples"))
