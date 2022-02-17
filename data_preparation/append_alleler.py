import os

##################
## Append Files ##
##################

# Creating a list of filenames
#filenames = ['file_1.txt', 'file_2.txt'] 
filenames = os.listdir('/Users/osman/Desktop/plot_test/alleler_files') 
# Open file3 in write mode
with open('file_3.txt', 'w') as outfile:
  
    # Iterate through list
    for names in filenames:
  
        # Open each file in read mode
        with open(names) as infile:
  
            # read the data from file1 and
            # file2 and write it in file3
            outfile.write(infile.read())
  
        # Add '\n' to enter data of file2
        # from next line
        outfile.write("")
