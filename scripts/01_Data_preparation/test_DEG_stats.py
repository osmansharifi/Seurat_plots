import os
import sys
import pandas as pd

def process_directory(directory_path):
    # Initialize a list to store the results
    results = []

    # Print the header for the results table
    print("Results:")
    print("Directory\tPositive_logFC\tNegative_logFC\tTotal_rows")

    # Iterate through each directory in the specified path
    for root, dirs, files in os.walk(directory_path):
        # Check if "DEGs.xlsx" file is in the current directory
        if "DEGs.xlsx" in files and "activated" not in os.path.basename(os.path.normpath(root)) and "not" not in os.path.basename(os.path.normpath(root)):
            file_path = os.path.join(root, "DEGs.xlsx")

            # Read the first 5 rows of the Excel file
            df = pd.read_excel(file_path, nrows=5)

            # Filter rows based on adj.P.Val <= 0.05
            significant_rows = df[df['adj.P.Val'] <= 0.05]

            # Count positive and negative logFC values
            positive_logFC_count = len(significant_rows[significant_rows['logFC'] > 0])
            negative_logFC_count = len(significant_rows[significant_rows['logFC'] < 0])

            # Total number of rows where adj.P.Val is <= 0.05
            total_rows = len(significant_rows)

            # Append results to the list
            results.append((os.path.basename(os.path.normpath(root)), positive_logFC_count, negative_logFC_count, total_rows))

            # Print the results for the current directory
            print(f"{os.path.basename(os.path.normpath(root))}\t{positive_logFC_count}\t{negative_logFC_count}\t{total_rows}")

    # Create a DataFrame from the results list
    df_results = pd.DataFrame(results, columns=['Directory', 'Positive_logFC', 'Negative_logFC', 'Total_rows'])

    # Save the DataFrame to a CSV file in the parent directory
    csv_output_path = os.path.join(directory_path, 'results.csv')
    df_results.to_csv(csv_output_path, index=False)
    print(f"\nResults saved to: {csv_output_path}")

if __name__ == "__main__":
    # Check if a command-line argument is provided
    if len(sys.argv) != 2:
        print("Usage: python script.py /path/to/parent/folder")
        sys.exit(1)

    # Get the parent directory path from the command-line argument
    parent_directory = sys.argv[1]

    # Call the function to process the directories and save results
    process_directory(parent_directory)
