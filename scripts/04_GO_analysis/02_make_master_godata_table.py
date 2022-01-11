#!/usr/bin/env python3

import pandas as pd

# Tissue types associated with each time point
tissue_types = ["WB", "CORT", "CORT", "CORT"]

# Cell types 
cell_types = ["L2_3_IT", "L6", "Sst", "L5", "L4", "Pvalb", "Sncg", "Non_neuronal", "Oligo", "Vip", "Lamp5", "Astro", "Peri", "Endo"]

# Gene ontology methods
go_onts = ["BP", "CC", "MF"]

# Make a master data frame containing all GO data for top 20, top 5, and top 3 GO terms
master_go_df_top20 = pd.DataFrame()
master_go_df_top5 = pd.DataFrame()
master_go_df_top3 = pd.DataFrame()

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
for time_point, meta, tissue in zip(time_points, meta_folders, tissue_types):
    for cell_type in cell_types:
        for ont in go_onts:
            # Read in GO Data
            data = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}/{cell_type}_M_MUT_and_WT_M_{time_point}_{tissue}_{ont}_gentable.csv')
            # Turn the GO Data into a pandas data frame
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

# Fill data frame with GO data (top 5)
for time_point, meta, tissue in zip(time_points, meta_folders, tissue_types):
    for cell_type in cell_types:
        for ont in go_onts:
            # Read in GO Data
            data = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}/{cell_type}_M_MUT_and_WT_M_{time_point}_{tissue}_{ont}_top5_gentable.csv')
            # Turn the GO Data into a pandas data frame
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
            master_go_df_top5 = master_go_df_top5.append(df)

# Fill data frame with GO data (top 3)
for time_point, meta, tissue in zip(time_points, meta_folders, tissue_types):
    for cell_type in cell_types:
        for ont in go_onts:
            # Read in GO Data
            data = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}/{cell_type}_M_MUT_and_WT_M_{time_point}_{tissue}_{ont}_top3_gentable.csv')
            # Turn the GO Data into a pandas data frame
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
            master_go_df_top3 = master_go_df_top3.append(df)

'''
#################
## For Females ##
#################

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
        for ont in go_onts:
            # Read in GO Data
            data = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}/{cell_type}_M_MUT_and_WT_F_{time_point}_{tissue}_{ont}_gentable.csv')
            # Turn the GO Data into a pandas data frame
            df = pd.DataFrame(data)
            
            # Add a new column to the data frame with the name of the data set
            metadata_name = (f'{cell_type}_F_{time_point}_{tissue}_{ont}')
            sex = []
            new_column_values = []
            cluster_names = []
            times = []
            tiss = []
            ontology = []
            for i in range(len(df["Term"])):
                new_column_values.append(metadata_name)
                sex.append("F")
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

# Fill data frame with GO data (top 5)
for time_point, meta, tissue in zip(time_points, meta_folders, tissue_types):
    for cell_type in cell_types:
        for ont in go_onts:
            # Read in GO Data
            data = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}/{cell_type}_M_MUT_and_WT_F_{time_point}_{tissue}_{ont}_top5_gentable.csv')
            # Turn the GO Data into a pandas data frame
            df = pd.DataFrame(data)
            
            # Add a new column to the data frame with the name of the data set
            metadata_name = (f'{cell_type}_F_{time_point}_{tissue}_{ont}')
            sex = []
            new_column_values = []
            cluster_names = []
            times = []
            tiss = []
            ontology = []
            for i in range(len(df["Term"])):
                new_column_values.append(metadata_name)
                sex.append("F")
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
            master_go_df_top5 = master_go_df_top5.append(df)

# Fill data frame with GO data (top 3)
for time_point, meta, tissue in zip(time_points, meta_folders, tissue_types):
    for cell_type in cell_types:
        for ont in go_onts:
            # Read in GO Data
            data = pd.read_csv(f'../../GO_data/GO_term_tables/{meta}/{cell_type}_M_MUT_and_WT_F_{time_point}_{tissue}_{ont}_top3_gentable.csv')
            # Turn the GO Data into a pandas data frame
            df = pd.DataFrame(data)
            
            # Add a new column to the data frame with the name of the data set
            metadata_name = (f'{cell_type}_F_{time_point}_{tissue}_{ont}')
            sex = []
            new_column_values = []
            cluster_names = []
            times = []
            tiss = []
            ontology = []
            for i in range(len(df["Term"])):
                new_column_values.append(metadata_name)
                sex.append("F")
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
            master_go_df_top3 = master_go_df_top3.append(df)
'''

###############################
## Create Master Data Frames ##
###############################

master_go_df_top20.to_csv("../../GO_data/GO_term_tables/master_go_data_males_top20.csv")
master_go_df_top5.to_csv("../../GO_data/GO_term_tables/master_go_data_males_top5.csv")   
master_go_df_top3.to_csv("../../GO_data/GO_term_tables/master_go_data_males_top3.csv")