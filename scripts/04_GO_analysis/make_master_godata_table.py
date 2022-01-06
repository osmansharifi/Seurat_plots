#!/usr/bin/env python3

import csv
import pandas as pd

# Folders containing input CSV data
meta_folders = ["M_MUT_and_WT_M_E18_WB",
            "M_MUT_and_WT_M_P30_CORT",
            "M_MUT_and_WT_M_P60_CORT",
            "M_MUT_and_WT_M_P120_CORT"]

# Time points
time_points = ["E18_WB", "P30_CORT", "P60_CORT", "P120_CORT"]

# Cell types 
cell_types = ["L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo"]

# Gene ontology methods
go_onts = ["BP", "CC", "MF"]

for tp in time_point:
    for cell_type in cell_types:
        for ont in go_onts:
        
        
        
for meta in metafolders: 
    for cell_type in zip(cell_types):
    
        df = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}




# Read in all GO Data from individual CSV files
with open

# Add a column containing the data name to each data frame
# Combine all data frames into a master data frame
# Output and save master data frame