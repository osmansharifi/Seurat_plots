import os
import pandas
import csv

# Path to file output
output_file_path = "."

# Mutation status
mut_stat = ["MUT", "WT"]

# Sex
sexes = ["M", "F"]

# Time Points
time_points = ["E18", "P30", "P60", "P120", "P150"]

# Tissues
tissue_types = ["WB", "CORT", "CORT", "CORT", "CORT"]

# Replicates
replicates = ["1", "2", "3", "4"]

# Make a master data frame to combine all the data from each individual sequence.alleler file
master_sequence_alleler = pd.DataFrame()

for mut in mut_stat:
    for sex in sexes:
        for time_point, tissue in zip(time_points, tissue_types):
            for rep in replicates:
                # Check if file exists (since males and females have different time points and time points may have different replicates)
                my_file = f'/share/lasallelab/Osman/cell_parsing_test/{mut}_{sex}_{time_point}_{tissue}{rep}/sequence.alleler'
                # Store True if file exists
                isExist = os.path.isfile(my_file)
                if isExist = 
                    with open(f'{my_file}') as fp:
                        # Read in the file contents
                        data = fp.read()
                        # Turn the file contents into a pandas data frame
                        df = pd.DataFrame(data)
                        # Append data to master data frame
                        master_sequence_alleler = master_sequence_alleler.append(df)
                else:
                    print(f'Skipping the following file because it cannot be found: {my_file}')

master_sequence_alleler.to_csv(f'{output_file_path}/master_sequence_alleler.csv')

# Convert CSV to TSV
with open(f'{output_file_path}/master_sequence_alleler.csv', 'r') as csvin, open(f'{output_file_path}/master_sequence_alleler.tsv') as tsvout:
    csvin = csv.reader(csvin)
    tsvout = csv.writer(tsvout, delimiter = '\t')
    
    for row in csvin:
        tsvout.writerow(row)