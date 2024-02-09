########################
## Upset plot of DEGs ##
########################
# Author : Osman Sharifi

## Load libraries
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

pdf(glue("{base_path}overlapping_maleMouse_human/allmouse_human_DEGs_upset.pdf"))
upset(fromList(listInput), sets = c('Male_mouse_Glut', 'Male_mouse_GABA','Female_mouse_Glut', 'Female_mouse_GABA', 'Female_human_Glut', 'Female_human_GABA'), keep.order = TRUE)
dev.off()


