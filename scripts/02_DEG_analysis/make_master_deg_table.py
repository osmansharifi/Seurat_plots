#!/usr/bin/env python3

import pandas as pd
import os
import sys

# Tissue types associated with each time point
tissue_types = ["WB", "CORT", "CORT", "CORT"]

# Cell types 
cell_types = ["L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo"]

# DEG methods
deg_tools = ["DEsingle", "LimmaVoomCC"]

# Make a master data frame containing all deg data for top 20 SYMBOLs
master_deg_df = pd.DataFrame()

###############
## For Males ##
###############

# Folders containing input CSV data
meta_folders = ["M_MUT_and_WT_M_E18_WB", "M_MUT_and_WT_M_P30_CORT", "M_MUT_and_WT_M_P60_CORT", "M_MUT_and_WT_M_P120_CORT"]

# Time points
time_points = ["E18", "P30", "P60", "P120"]

# Fill data frame with deg data (top 20)
for meta, deg, cell_type, time_point, tissue_type in zip(meta_folders, deg_tools, cell_types, time_points, tissue_types):
	for meta in meta_folders:
		for deg in deg_tools:
			for cell_type in cell_types:
				for time_point in time_points:
					my_path = f'/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/{meta}/{deg}/{cell_type}/DEGs.xlsx'
					if not os.path.isfile(my_path):
						print("file was not found", my_path)	 
						continue
					print(f'Adding data from {my_path} to data frame')
					# Read in DEG Data
					data = pd.read_excel(my_path)
					# Turn the DEG Data into a pandas data frame
					df = pd.DataFrame(data)
			
					# Add a new column to the data frame with the name of the data set
					metadata_name = (f'{meta}')
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
						tiss.append(f'{tissue_type}')
						deg_method.append(f'{deg}')
					df.insert(1, "Metadata", new_column_values, True)
					df.insert(2, "Sex", sex, True)
					df.insert(3, "Cell Type", cluster_names, True)
					df.insert(4, "Time Point", times, True)
					df.insert(5, "Tissue", tiss, True)
					df.insert(6, "deg_method", deg_method, True)
					master_deg_df = master_deg_df.append(df)

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
for meta, deg, cell_type, time_point, tissue_type in zip(meta_folders, deg_tools, cell_types, time_points, tissue_types):
	for meta in meta_folders:
		for deg in deg_tools:
			for cell_type in cell_types:
				my_path = f'/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/{meta}/{deg}/{cell_type}/DEGs.xlsx'
				if not os.path.isfile(my_path):
					print("ERROR: file was not found", my_path)	 
					continue
				print(f'Adding data from {my_path} to data frame')
				# Read in DEG Data
				data = pd.read_excel(my_path)
				# Turn the DEG Data into a pandas data frame
				df = pd.DataFrame(data)
			
				# Add a new column to the data frame with the name of the data set
				metadata_name = (f'{meta}')
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
					tiss.append(f'{tissue_type}')
					deg_method.append(f'{deg}')
				df.insert(1, "Metadata", new_column_values, True)
				df.insert(2, "Sex", sex, True)
				df.insert(3, "Cell Type", cluster_names, True)
				df.insert(4, "Time Point", times, True)
				df.insert(5, "Tissue", tiss, True)
				df.insert(6, "deg_method", deg_method, True)
				master_deg_df = master_deg_df.append(df)
             
###############################
## Create Master Data Frames ##
###############################

master_deg_df.to_csv("/Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/master_deg_data_top20.csv")
