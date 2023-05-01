#################################
## Create a df from xlsx files ##
#################################

library(dplyr)

#Load excel file paths
dir_path_30_f <- "/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/females/M_MUT_and_WT_F_P30_CORT/limmaVoomCC"
dir_path_60_f <- "/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/females/M_MUT_and_WT_F_P60_CORT/limmaVoomCC"
dir_path_150_f <- "/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/females/M_MUT_and_WT_F_P150_CORT/limmaVoomCC"
dir_path_30_m <- "/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/males/M_MUT_and_WT_M_P30_CORT/limmaVoomCC"
dir_path_60_m <- "/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/males/M_MUT_and_WT_M_P60_CORT/limmaVoomCC"
dir_path_120_m <- "/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/males/M_MUT_and_WT_M_P120_CORT/limmaVoomCC"
dir_path_human <- "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/H_MUT_and_WT_F_ALL_CORT/limmaVoomCC"
dir_paths = c(dir_path_human) # dir_path_30_m, dir_path_60_m, dir_path_120_m, dir_path_30_f, dir_path_60_f, dir_path_150_f

# Create an empty list to store the file names and contents
file_contents <- list()

# Loop through each directory and read the Excel file
for (dir_path in dir_paths) {
  # Get the list of directories
  dir_list <- list.dirs(dir_path, recursive = FALSE)
  dir_list <- dir_list[-c(1, which(grepl("interactivePlots", dir_list)=="TRUE"), which(grepl("-activated", dir_list)=="TRUE"), which(grepl("-not", dir_list)=="TRUE"))]
  
  # Loop through each directory and read the Excel file
  for (dir_name in dir_list) {
    # Get the file path
    file_path <- file.path(dir_name, "enrichr.xlsx")
    
    # Read the Excel file if it exists
    if (file.exists(file_path)) {
      # Read the Excel file and store the contents in a data frame
      file_contents[[basename(dir_name)]] <- readxl::read_excel(file_path, sheet = "KEGG_2019_Mouse")
    }
  }
  
  # Create a column in the dataframe and add Cell_Type info
  for (i in seq_along(file_contents)) {
    file_contents[[i]]$Cell_Type <- names(file_contents)[i]
  }
  
  #Create a dataframe with all the KEGG terms
  mouse <- bind_rows(file_contents, .id = "Cell_Type")
  
  # Extract the part of the string between "Differential_expression/" and "M_MUT_and_WT_F_P30_CORT/"
  mystring <- gsub(".*\\/Differential_expression\\/", "", dir_path)
  mystring <- gsub("\\/limmaVoomCC.*", "", mystring)
  
  # Split each string into two parts
  mystrings_split <- strsplit(mystring, "/")
  
  # Add two new columns to the existing data frame
  mouse$Sex <- mystrings_split[[1]][1]
  mouse$metadata <- mystrings_split[[1]][2]
  
  # Append to previous data
  if(exists("full_mouse")) {
    full_mouse <- rbind(full_mouse, mouse)
  } else {
    full_mouse <- mouse
  }
  
  # Clear the list for the next iteration
  file_contents <- list()
}

# assign full_mouse to a variable that represent the loaded samples 
male_30 <- full_mouse

# Concatenate the three data frames
concatenated_df <- bind_rows(full_female, male_30, male_60_120)

#save dataframe as csv
write.csv(concatenated_df, file = "mouse_total_kegg.csv", row.names = FALSE)


filtered_mouse <- concatenated_df[concatenated_df$Adjusted.P.value <= 0.05,]
# Write the merged_data dataframe to a CSV file
write.csv(filtered_mouse, file = "mouse_total_sig_kegg.csv", row.names = FALSE)

#human kegg terms
human_kegg <- full_mouse
write.csv(human_kegg, file = "human_total_kegg.csv", row.names = FALSE)
filtered_human <- human_kegg[human_kegg$Adjusted.P.value <= 0.05,]
# Write the merged_data dataframe to a CSV file
write.csv(filtered_human, file = "human_total_sig_kegg.csv", row.names = FALSE)
