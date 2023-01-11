# This is code to concatnate single cell DEGs from all cell types and create a csv containing all DEG info and metadata

# Load libraries
library(dplyr)
library(tidyverse)
library(glue)
library(readxl)

# Cell types 
cell_types = c("L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo")

# DEG methods
deg_tools = c("LimmaVoomCC")

# Path to files
base_path = ("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression")
# Create a vector containing all the files 
e18_males <- list.files(glue("{base_path}/males/M_MUT_and_WT_M_E18_WB/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
p30_males <- list.files(glue("{base_path}/males/M_MUT_and_WT_M_P30_CORT/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
p60_males <- list.files(glue("{base_path}/males/M_MUT_and_WT_M_P60_CORT/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
p120_males <- list.files(glue("{base_path}/males/M_MUT_and_WT_M_P120_CORT/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
e18_females <- list.files(glue("{base_path}/females/M_MUT_and_WT_F_E18_WB/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
p30_females <- list.files(glue("{base_path}/females/M_MUT_and_WT_F_P30_CORT/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
p60_females <- list.files(glue("{base_path}/females/M_MUT_and_WT_F_P60_CORT/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
p150_females <- list.files(glue("{base_path}/females/M_MUT_and_WT_F_P150_CORT/{deg_tools}/{cell_types}"), pattern = "DEGs.xlsx", full.names = TRUE)
all_files_paths <- c(e18_males, p30_males, p60_males, p120_males, e18_females, p30_females, p60_females, p150_females)

# List of file locations
file_locations <- all_files_paths[-grep("sig_DEGs.xlsx", all_files_paths, fixed=T)]

# Create an empty dataframe
df_final <- data.frame()

# Iterate over each file location
for(i in 1:length(file_locations)){
  # Read the Excel file into a dataframe
  file_path <- paste0("",file_locations[i])
  df <- read_excel(file_path)
  
  # Add a new column for the file path
  df$file_path <- file_path
  
  # Append the dataframe to the final dataframe
  df_final <- rbind(df_final, df)
}

# Remove unimportant columns
df_final$t <- NULL
df_final$B <- NULL

# Add metadata columns
split_file_path <- str_split_fixed(df_final$file_path, "/", n = 12)
df_final$Sex <- split_file_path[,8]
df_final$DEG_method <- split_file_path[,10]
df_final$Cell_Type <- split_file_path[,11]
df_final$Metadata <- split_file_path[,9]
df_final$file_path <- NULL

split_metadata <- str_split_fixed(df_final$Metadata, "_", n = 7)
df_final$Time_point <- split_metadata[,6]
df_final$Tissue <- split_metadata[,7]
print(df_final)

# Write the final dataframe to a new csv file
write.csv(df_final,"{base_path}/total_mouse_DEGs_limmaVoom.csv",row.names = TRUE)

# Write csv containing only the significant DEGs based on adj.P.Val
total_sig_mouse_DEGs_limmaVoom <- filter(df_final, adj.P.Val <= 0.05)
write.csv(total_sig_mouse_DEGs_limmaVoom,"{base_path}/total_sig_mouse_DEGs_limmaVoom.csv", row.names = TRUE)
