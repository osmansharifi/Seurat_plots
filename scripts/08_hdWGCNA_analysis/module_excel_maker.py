import pandas as pd

# Load the CSV file
csv_file_path = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/enrich_df.csv"
df = pd.read_csv(csv_file_path)

# Filter rows where 'db' column contains 'GO_Biological_Process_2023'
filtered_df = df[df['db'].str.contains('KEGG_2019_Mouse', case=False, na=False)]

# Save the filtered DataFrame to a new CSV file
filtered_csv_file_path = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/filtered_enrich_df.csv"
filtered_df.to_csv(filtered_csv_file_path, index=False)

print(f"Filtered DataFrame saved to: {filtered_csv_file_path}")
