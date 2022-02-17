#!/usr/bin/env python3

import pandas as pd
import os

# Tissue types associated with each time point
file_names = ["file_1.txt", "file_2.txt"]

# Make a master data frame containing all GO data for top 20, top 5, and top 3 GO terms
master_mecp2_df = pd.DataFrame()

for lines in (file_names):
	my_path = '/Users/osman/Desktop/plot_test/alleler_files/file_3.txt'
		data = pd.read_csv(my_path, sep='')
		print(data)

'''
    df = pd.DataFrame(data)
                
                # Add a new column to the data frame with the name of the data set
                metadata_name = (f'{cell_type}_M_{time_point}_{tissue}_{ont}')
                sex = []
                new_column_values = []
                cluster_names = []
                times = []
                tiss = []
                ontology = []
                for i in range(len(df["Term"])):
                    new_column_values.append(metadata_name)
                    sex.append("M")
                    cluster_names.append(f'{cell_type}')
                    times.append(f'{time_point}')
                    tiss.append(f'{tissue}')
                    ontology.append(f'{ont}')
                df.insert(1, "Metadata", new_column_values, True)
                df.insert(2, "Sex", sex, True)
                df.insert(3, "Cell Type", cluster_names, True)
                df.insert(4, "Time Point", times, True)
                df.insert(5, "Tissue", tiss, True)
                df.insert(6, "Ontology", ontology, True)
                master_go_df_top20 = master_go_df_top20.append(df)
            else:
                print(f'Skipping the following file because it cannot be found: {my_path}')
'''