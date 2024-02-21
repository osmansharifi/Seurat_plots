import sys
import pandas as pd

def analyze_csv_file(file_path):
    try:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(file_path)

        # Filter rows based on the condition in the adj.P.Val column
        filtered_df = df[df['adj.P.Val'] <= 0.05]

        # Count positive and negative values in the logFC column
        positive_count = (filtered_df['logFC'] > 0).sum()
        negative_count = (filtered_df['logFC'] < 0).sum()

        # Display the results
        print(f"Number of rows with adj.P.Val <= 0.05: {len(filtered_df)}")
        print(f"Number of rows with positive logFC values: {positive_count}")
        print(f"Number of rows with negative logFC values: {negative_count}")

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
    except pd.errors.EmptyDataError:
        print(f"Error: The CSV file at {file_path} is empty or not properly formatted")

if __name__ == "__main__":
    # Check if the file path is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python script.py /path/to/your/file.csv")
    else:
        # Get the file path from the command-line argument
        file_path = sys.argv[1]

        # Call the function to analyze the CSV file
        analyze_csv_file(file_path)
