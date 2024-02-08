########################
## Upset plot of DEGs ##
########################
# Author : Osman Sharifi

# Load libraries
library(UpSetR)
library(dplyr)
library(ggplot2)
library(glue)

#####################
## Load mouse DEGs ##
#####################
testing <- read.csv('/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/total_sig_mouse_DEGs_limmaVoom.csv')
mouse_DEGs <- read.csv('/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/total_mouse_DEGs_limmaVoom.csv') %>%
  select(-X) %>%
  mutate(
    broad_class = case_when(
      Cell_Type %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Vip") ~ "GABAergic",
      Cell_Type %in% c("L2_3_IT", "L4", "L5", "L6") ~ "Glutamatergic",
      Cell_Type %in% c("Astro", "Non-neuronal", "Oligo") ~ "Non-neuronal",
      TRUE ~ "Other"
    )
  )

mouse_male_DEGs <- mouse_DEGs %>%
  filter(Sex != 'females' & Time_point != 'E18' & adj.P.Val <= 0.05) %>%
  mutate(SYMBOL = toupper(SYMBOL))

mouse_female_DEGs <- mouse_DEGs %>%
  filter(Sex != 'males' & Time_point != 'E18' & adj.P.Val <= 0.05) %>%
  mutate(SYMBOL = toupper(SYMBOL))

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
#Glut_human_female_DEGs <- human_female_DEGs$SYMBOL[human_female_DEGs$broad_class == 'Glutamatergic']
#Glut_mouse_male_DEGs <- mouse_male_DEGs$SYMBOL[mouse_male_DEGs$broad_class == 'Glutamatergic']
#Glut_mouse_female_DEGs <- mouse_male_DEGs$SYMBOL[mouse_male_DEGs$broad_class == 'Glutamatergic']
#GABA_human_female_DEGs <- human_female_DEGs$SYMBOL[human_female_DEGs$broad_class == 'GABAergic']
#GABA_mouse_male_DEGs <- mouse_male_DEGs$SYMBOL[mouse_male_DEGs$broad_class == 'GABAergic']
#GABA_mouse_female_DEGs <- mouse_female_DEGs$SYMBOL[mouse_female_DEGs$broad_class == 'GABAergic']
#Non_neuronal_human_female_DEGs <- human_female_DEGs$SYMBOL[human_female_DEGs$broad_class == 'Non-neuronal']
#Non_neuronal_mouse_male_DEGs <- mouse_male_DEGs$SYMBOL[mouse_male_DEGs$broad_class == 'Non-neuronal']
#Non_neuronal_mouse_female_DEGs <- mouse_female_DEGs$SYMBOL[mouse_female_DEGs$broad_class == 'Non-neuronal']
#intersection_Glut <- intersect(Glut_human_female_DEGs,Glut_mouse_male_DEGs)
#intersection_GABA <- intersect(GABA_human_female_DEGs,GABA_mouse_male_DEGs)
#intersection_non_neuronal <- intersect(Non_neuronal_human_female_DEGs,Non_neuronal_mouse_male_DEGs)
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/'

#######################
## Create upset plot ##
#######################
# Function to create a list from Time_Point and SYMBOL columns
create_list_from_df <- function(df) {
  time_points <- unique(df$Time_point)
  result_list <- lapply(time_points, function(tp) {
    SYMBOLs <- df$SYMBOL[df$Time_point == tp]
    return(setNames(SYMBOLs, NULL))  # Exclude "Time_Point" and "SYMBOL" as names
  })
  names(result_list) <- time_points
  return(result_list)
}

# Create lists for each dataframe
human_female_glut_list <- create_list_from_df(human_female_DEGs %>% filter(broad_class == "Glutamatergic"))
human_female_gaba_list <- create_list_from_df(human_female_DEGs %>% filter(broad_class == "GABAergic"))
mouse_female_glut_list <- create_list_from_df(mouse_female_DEGs %>% filter(broad_class == "Glutamatergic"))
mouse_female_gaba_list <- create_list_from_df(mouse_female_DEGs %>% filter(broad_class == "GABAergic"))
mouse_male_glut_list <- create_list_from_df(mouse_male_DEGs %>% filter(broad_class == "Glutamatergic"))
mouse_male_gaba_list <- create_list_from_df(mouse_male_DEGs %>% filter(broad_class == "GABAergic"))

#Create input lists and make upset PDF
listInput <- list(Female_human_Glut = human_female_glut_list$Postmortem, 
                  Female_human_GABA = human_female_gaba_list$Postmortem, 
                  Female_mouse_Glut = mouse_female_glut_list$P30,
                  Female_mouse_GABA = mouse_female_gaba_list$P30, 
                  Male_mouse_Glut = mouse_male_glut_list$P30, 
                  Male_mouse_GABA = mouse_male_gaba_list$P30)
pdf(glue("{base_path}overlapping_maleMouse_human/Glut_allmouse_human_upset.pdf"))
upset(fromList(listInput), sets = c('Female_human_Glut', 'Female_mouse_Glut',' Male_mouse_Glut'), keep.order = TRUE)
dev.off()


listInput <- list(Glut_human_female_DEGs, 
                  GABA_human_female_DEGs, 
                  Glut_mouse_female_DEGs,
                  GABA_mouse_female_DEGs, 
                  Glut_mouse_male_DEGs, 
                  Glut_mouse_male_DEGs)
pdf(glue("{base_path}overlapping_maleMouse_human/GABA_allmouse_human_upset.pdf"))
upset(fromList(listInput), sets = c('Female_human_Glut', 'Female_mouse_Glut',' Male_mouse_Glut'), keep.order = TRUE)
dev.off()

upset(fromList(listInput),keep.order = TRUE)
################################################
# Load libraries
library(UpSetR)
library(dplyr)
library(ggplot2)
library(glue)
library(readxl)  # Added for read_xlsx function

# Function to read and filter KEGG terms
read_and_filter_kegg <- function(file_path) {
  kegg_terms <- readxl::read_xlsx(file_path, col_names = TRUE)
  kegg_terms <- filter(kegg_terms, P.value <= 0.05)
  
  # Convert column name to "broad_class" (case-insensitive)
  col_name <- tolower("broad_class")
  if (col_name %in% names(kegg_terms)) {
    kegg_terms$broad_class <- kegg_terms[[col_name]]
    kegg_terms <- select(kegg_terms, -col_name)
  }
  
  return(kegg_terms)
}

# Function to create a list from broad_class and SYMBOL columns
create_list_from_df <- function(df) {
  broad_classes <- unique(df$broad_class)
  result_list <- lapply(broad_classes, function(bc) {
    symbols <- df$SYMBOL[df$broad_class == bc]
    return(setNames(symbols, NULL))  # Exclude "broad_class" and "SYMBOL" as names
  })
  names(result_list) <- broad_classes
  return(result_list)
}

# Set base path
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/'

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

# Create lists for each dataframe
expr3_gaba_list <- create_list_from_df(expr3_gaba)
expr3_glut_list <- create_list_from_df(expr3_glut)
expr4_gaba_list <- create_list_from_df(expr4_gaba)
expr4_glut_list <- create_list_from_df(expr4_glut)

# Create input lists and make upset PDF for Glut
listInput_glut <- list(Female_human = expr3_glut_list, 
                       Female_mouse = expr4_glut_list, 
                       Male_mouse = expr4_glut_list)

pdf(glue("{base_path}overlapping_maleMouse_human/Glut_allmouse_human_upset.pdf"))
upset(fromList(listInput_glut), sets = c('Female_human', 'Female_mouse', 'Male_mouse'), keep.order = TRUE)
dev.off()

# Create input lists and make upset PDF for GABA
listInput_gaba <- list(Female_human = expr3_gaba_list, 
                       Female_mouse = expr4_gaba_list, 
                       Male_mouse = expr4_gaba_list)

pdf(glue("{base_path}overlapping_maleMouse_human/GABA_allmouse_human_upset.pdf"))
upset(fromList(listInput_gaba), sets = c('Female_human', 'Female_mouse', 'Male_mouse'), keep.order = TRUE)
dev.off()



################



# Load libraries
library(UpSetR)
library(dplyr)
library(glue)

# Load mouse DEGs
mouse_DEGs <- read.csv('/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/total_mouse_DEGs_limmaVoom.csv') %>%
  select(-X) %>%
  mutate(
    broad_class = case_when(
      Cell_Type %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Vip") ~ "GABAergic",
      Cell_Type %in% c("L2_3_IT", "L4", "L5", "L6") ~ "Glutamatergic",
      Cell_Type %in% c("Astro", "Non-neuronal", "Oligo") ~ "Non-neuronal",
      TRUE ~ "Other"
    )
  )

mouse_male_DEGs <- mouse_DEGs %>%
  filter(Sex != 'females' & Time_point != 'E18' & adj.P.Val <= 0.05) %>%
  mutate(SYMBOL = toupper(SYMBOL))

mouse_female_DEGs <- mouse_DEGs %>%
  filter(Sex != 'males' & Time_point != 'E18' & adj.P.Val <= 0.05) %>%
  mutate(SYMBOL = toupper(SYMBOL))

# Load human DEGs
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

# Create upset plot
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/'

# Function to create a list from Time_Point and SYMBOL columns
create_list_from_df <- function(df) {
  time_points <- unique(df$Time_point)
  result_list <- lapply(time_points, function(tp) {
    SYMBOLs <- df$SYMBOL[df$Time_point == tp]
    return(setNames(SYMBOLs, NULL))  # Exclude "Time_Point" and "SYMBOL" as names
  })
  names(result_list) <- time_points
  return(result_list)
}

# Create lists for each dataframe
human_female_glut_list <- create_list_from_df(human_female_DEGs %>% filter(broad_class == "Glutamatergic"))
human_female_gaba_list <- create_list_from_df(human_female_DEGs %>% filter(broad_class == "GABAergic"))
mouse_female_glut_list <- create_list_from_df(mouse_female_DEGs %>% filter(broad_class == "Glutamatergic"))
mouse_female_gaba_list <- create_list_from_df(mouse_female_DEGs %>% filter(broad_class == "GABAergic"))
mouse_male_glut_list <- create_list_from_df(mouse_male_DEGs %>% filter(broad_class == "Glutamatergic"))
mouse_male_gaba_list <- create_list_from_df(mouse_male_DEGs %>% filter(broad_class == "GABAergic"))

# Create input lists and make upset PDF
listInput <- list(Female_human_Glut = human_female_glut_list$Postmortem, 
                  Female_human_GABA = human_female_gaba_list$Postmortem, 
                  Female_mouse_Glut = c(mouse_female_glut_list$P30, mouse_female_glut_list$P60),
                  Female_mouse_GABA = c(mouse_female_gaba_list$P30, mouse_female_gaba_list$P150),
                  Male_mouse_Glut = c(mouse_male_glut_list$P30, mouse_male_glut_list$P60, mouse_male_glut_list$P120),
                  Male_mouse_GABA = c(mouse_male_gaba_list$P30, mouse_male_gaba_list$P60, mouse_male_gaba_list$P120))

pdf(glue("{base_path}overlapping_maleMouse_human/Glut_allmouse_human_upset.pdf"))
upset(fromList(listInput), sets = c('Female_human_Glut', 'Female_human_GABA', 'Female_mouse_Glut', 'Female_mouse_GABA', 'Male_mouse_Glut', 'Male_mouse_GABA'), keep.order = TRUE)
dev.off()
