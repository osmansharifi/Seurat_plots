########################
## Upset plot of DEGs ##
########################
# Author : Osman Sharifi

# Load libraries
library(UpSetR)
library(dplyr)
library(ggplot2)
library(glue)

# Load data
read_and_filter_kegg <- function(file_path) {
  kegg_terms <- readxl::read_xlsx(file_path, col_names = TRUE)
  kegg_terms <- filter(kegg_terms, P.value <= 0.05)
  
  # Convert column name to "Time_Point" (case-insensitive)
  col_name <- tolower("Time_Point")
  if (col_name %in% names(kegg_terms)) {
    kegg_terms$Time_Point <- kegg_terms[[col_name]]
    kegg_terms <- select(kegg_terms, -col_name)
  }
  
  return(kegg_terms)
}
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'

# Read and filter GABA kegg terms
gaba_kegg_terms_path <- glue::glue("{base_path}broad_group_analysis/GABA_kegg_all.xlsx")
expr3_gaba <- read_and_filter_kegg(gaba_kegg_terms_path)

# Read and filter Glut kegg terms
glut_kegg_terms_path <- glue::glue("{base_path}broad_group_analysis/Glut_kegg_all.xlsx")
expr3_glut <- read_and_filter_kegg(glut_kegg_terms_path)

# Read and filter mut_from_mut_vs_wt_from_wt Glut kegg terms
mut_glut_kegg_terms_path <- glue::glue("{base_path}mut_from_mut_vs_wt_from_wt/Glut_kegg_all.xlsx")
expr4_glut <- read_and_filter_kegg(mut_glut_kegg_terms_path)

# Read and filter mut_from_mut_vs_wt_from_wt GABA kegg terms
mut_gaba_kegg_terms_path <- glue::glue("{base_path}mut_from_mut_vs_wt_from_wt/GABA_kegg_all.xlsx")
expr4_gaba <- read_and_filter_kegg(mut_gaba_kegg_terms_path)

#######################
## Create upset plot ##
#######################
# Function to create a list from Time_Point and Term columns
create_list_from_df <- function(df) {
  time_points <- unique(df$Time_point)
  result_list <- lapply(time_points, function(tp) {
    terms <- df$Term[df$Time_point == tp]
    return(setNames(terms, NULL))  # Exclude "Time_Point" and "Term" as names
  })
  names(result_list) <- time_points
  return(result_list)
}

# Create lists for each dataframe
expr3_gaba_list <- create_list_from_df(expr3_gaba)
expr3_glut_list <- create_list_from_df(expr3_glut)
expr4_gaba_list <- create_list_from_df(expr4_gaba)
expr4_glut_list <- create_list_from_df(expr4_glut)

#Create input lists and make upset PDF
listInput <- list(Experiment3_P30 = expr3_glut_list$P30, 
                  Experiment3_P60 = expr3_glut_list$P60, 
                  Experiment3_P150 = expr3_glut_list$P150,
                  Experiment4_P30 = expr4_glut_list$P30, 
                  Experiment4_P60 = expr4_glut_list$P60, 
                  Experiment4_P150 = expr4_glut_list$P150)
pdf(glue("{base_path}upset_glut_kegg.pdf"))
upset(fromList(listInput), sets = c('Experiment4_P150','Experiment3_P150', 'Experiment4_P60','Experiment3_P60', 'Experiment4_P30', 'Experiment3_P30'), keep.order = TRUE)
dev.off()

listInput <- list(Experiment3_P30 = expr3_gaba_list$P30,
                  Experiment3_P60 = expr3_gaba_list$P60,
                  Experiment3_P150 = expr3_gaba_list$P150,
                  Experiment4_P30 = expr4_gaba_list$P30,
                  Experiment4_P60 = expr4_gaba_list$P60,
                  Experiment4_P150 = expr4_gaba_list$P150)
pdf(glue("{base_path}upset_gaba_kegg.pdf"))
upset(fromList(listInput), sets = c('Experiment4_P150','Experiment3_P150', 'Experiment4_P60','Experiment3_P60', 'Experiment4_P30', 'Experiment3_P30'), keep.order = TRUE)
dev.off()
