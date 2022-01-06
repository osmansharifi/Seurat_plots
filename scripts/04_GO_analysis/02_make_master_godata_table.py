#!/usr/bin/env python3

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

# Make a master data frame containing all GO data
master_go_df = pd.DataFrame()

# Fill data frame with GO data
for time_point, meta in zip(time_points, meta_folders):
    for cell_type in cell_types:
        for ont in go_onts:
            # Read in GO Data
            data = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}/{cell_type}_M_MUT_and_WT_M_{time_point}_{ont}_gentable.csv')
            # Turn the GO Data into a pandas data frame
            df = pd.DataFrame(data)
            
            # Add a new column to the data frame with the name of the data set
            metadata_name = (f'{cell_type}_{time_point}_{ont}')
            new_column_values = []
            for i in range(len(df["Term"])):
                new_column_values.append(metadata_name)
            df.insert(0, "Metadata", new_column_values, True)
            master_go_df = master_go_df.append(df)

master_go_df.to_csv("../../GO_data/GO_term_tables/master_go_data_males.csv")      