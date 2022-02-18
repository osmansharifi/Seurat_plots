import os
import csv

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

# Make a list containing path names for files that exist
my_files = []

for mut in mut_stat:
    for sex in sexes:
        for time_point, tissue in zip(time_points, tissue_types):
            for rep in replicates:
                # Check if file exists (since males and females have different time points and time points may have different replicates)
                my_file = f'/share/lasallelab/Osman/cell_parsing_test/{mut}_{sex}_{time_point}_{tissue}{rep}/sequence.alleler'
                # Store True if file exists
                isExist = os.path.isfile(my_file)
                if isExist:
                    my_files.append(my_file)
                else:
                    print(f'Skipping the following file because it cannot be found: {my_file}')

with open('master_sequence_alleler.txt', 'w') as outfile:
    # Iterate through list
    for file in my_files:
        # Open each file in read mode
        with open(file) as infile:
            outfile.write(infile.read())

# Close file to avoid ValueError: I/O operation
outfile.close()