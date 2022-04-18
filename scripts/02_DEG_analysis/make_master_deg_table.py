#!/usr/bin/env python3

import pandas as pd
import os

# Tissue types associated with each time point
tissue_types = ["WB", "CORT", "CORT", "CORT"]

# Cell types 
cell_types = ["L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo"]

# DEG methods
deg_tool = ["DEsingle", "LimmaVoomCC"]

# Make a master data frame containing all GO data for top 20, top 5, and top 3 GO SYMBOLs
master_deg_df_top20 = pd.DataFrame()
master_deg_df_top5 = pd.DataFrame()
master_deg_df_top3 = pd.DataFrame()


###############
## For Males ##
###############

# Folders containing input CSV data
meta_folders = ["M_MUT_and_WT_M_E18_WB",
            "M_MUT_and_WT_M_P30_CORT",
            "M_MUT_and_WT_M_P60_CORT",
            "M_MUT_and_WT_M_P120_CORT"]

# Time points
time_points = ["E18", "P30", "P60", "P120"]

# Fill data frame with GO data (top 20)
for time_point, meta, tissue, deg in zip(time_points, meta_folders, tissue_types, deg_tool):
    for cell_type in cell_types:
        for deg in deg_tool:
            my_path = f'/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/{meta}/{deg_tool}/{cell_types}/DEGs.xlsx'
            isExist = os.path.isfile(my_path)
            if isExist:
                print(f'Adding data from {my_path} to data frame')
                # Read in GO Data
                data = pd.read_csv(my_path)
                # Turn the GO Data into a pandas data frame
                df = pd.DataFrame(data)
                
                # Add a new column to the data frame with the name of the data set
                metadata_name = (f'{cell_type}_M_{time_point}_{tissue}_{ont}')
                sex = []
                new_column_values = []
                cluster_names = []
                times = []
                tiss = []
                deg_method = []
                for i in range(len(df["SYMBOL"])):
                    new_column_values.append(metadata_name)
                    sex.append("M")
                    cluster_names.append(f'{cell_type}')
                    times.append(f'{time_point}')
                    tiss.append(f'{tissue}')
                    deg_method.append(f'{ont}')
                df.insert(1, "Metadata", new_column_values, True)
                df.insert(2, "Sex", sex, True)
                df.insert(3, "Cell Type", cluster_names, True)
                df.insert(4, "Time Point", times, True)
                df.insert(5, "Tissue", tiss, True)
                df.insert(6, "deg_method", deg_method, True)
                master_deg_df_top20 = master_deg_df_top20.append(df)
            else:
                print(f'Skipping the following file because it cannot be found: {my_path}')


#################
## For Females ##
#################

# Note: Some DEG data files are missing for females, so missing data are skipped over

# Folders containing input CSV data
meta_folders = ["M_MUT_and_WT_F_E18_WB",
            "M_MUT_and_WT_F_P30_CORT",
            "M_MUT_and_WT_F_P60_CORT",
            "M_MUT_and_WT_F_P150_CORT"]

# Time points
time_points = ["E18", "P30", "P60", "P150"]

# Fill data frame with GO data (top 20)
for time_point, meta, tissue in zip(time_points, meta_folders, tissue_types):
    for cell_type in cell_types:
        for ont in deg_tool:
            my_path = f'/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/{meta}/{deg_tool}/{cell_types}/DEGs.xlsx'
            isExist = os.path.isfile(my_path)
            if isExist:
                print(f'Adding data from {my_path} to data frame')
                # Read in GO Data
                data = pd.read_csv(my_path)
                # Turn the GO Data into a pandas data frame
                df = pd.DataFrame(data)
                
                # Add a new column to the data frame with the name of the data set
                metadata_name = (f'{cell_type}_F_{time_point}_{tissue}_{ont}')
                sex = []
                new_column_values = []
                cluster_names = []
                times = []
                tiss = []
                deg_method = []
                for i in range(len(df["SYMBOL"])):
                    new_column_values.append(metadata_name)
                    sex.append("F")
                    cluster_names.append(f'{cell_type}')
                    times.append(f'{time_point}')
                    tiss.append(f'{tissue}')
                    deg_method.append(f'{ont}')
                df.insert(1, "Metadata", new_column_values, True)
                df.insert(2, "Sex", sex, True)
                df.insert(3, "Cell Type", cluster_names, True)
                df.insert(4, "Time Point", times, True)
                df.insert(5, "Tissue", tiss, True)
                df.insert(6, "deg_method", deg_method, True)
                master_deg_df_top20 = master_deg_df_top20.append(df)
            else:
                print(f'Skipping the following file because it cannot be found: {my_path}')
                
###############################
## Create Master Data Frames ##
###############################

master_deg_df_top20.to_csv("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/master_deg_data_top20.csv")
