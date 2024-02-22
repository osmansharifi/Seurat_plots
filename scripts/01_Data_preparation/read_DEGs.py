import pandas as pd
import os
import sys

def count_symbols(df):
    # Filter rows where adj.P.Val <= 0.05
    significant_rows = df[df['adj.P.Val'] <= 0.05]

    # Count total number of SYMBOL
    total_symbols = len(significant_rows['SYMBOL'])

    # Count number of SYMBOL where logFC is positive
    positive_symbols = len(significant_rows[significant_rows['logFC'] > 0]['SYMBOL'])

    # Count number of SYMBOL where logFC is negative
    negative_symbols = len(significant_rows[significant_rows['logFC'] < 0]['SYMBOL'])

    return total_symbols, positive_symbols, negative_symbols

def process_directory(parent_directory, file_name):
    # Read the CSV file into a DataFrame
    file_path = os.path.join(parent_directory, file_name)
    df = pd.read_csv(file_path)

    # Get unique values in the 'Time_point' column
    time_points = df['Time_point'].unique()

    # Create a DataFrame to store the results
    df_results = pd.DataFrame(columns=['Time_point', 'Positive_logFC', 'Negative_logFC', 'Total_rows'])

    for time_point in time_points:
        # Filter DataFrame for the current time_point
        time_point_df = df[df['Time_point'] == time_point]

        # Get counts for the current time_point
        total_symbols, positive_symbols, negative_symbols = count_symbols(time_point_df)

        # Append results to the DataFrame
        df_results = df_results.append({
            'Time_point': time_point,
            'Positive_logFC': positive_symbols,
            'Negative_logFC': negative_symbols,
            'Total_rows': total_symbols
        }, ignore_index=True)

    return df_results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python read_DEGs.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]

    # Extract the directory and file name
    parent_directory, file_name = os.path.split(file_path)

    # Process the file and get the results
    df_results = process_directory(parent_directory, file_name)

    # Print the results
    print("Results:")
    print(df_results)

    # Save the results to a CSV file in the same directory as the input file
    csv_output_path = os.path.join(parent_directory, 'results.csv')
    df_results.to_csv(csv_output_path, index=False)
    print(f"Results saved to {csv_output_path}")
    