import pandas as pd

# Load the file
df = pd.read_csv("/Users/samuelgair/Desktop/BA_code/snakemake_fun/BAYESPRISM-updated_sig_matrix1002.txt", sep="\t", index_col=0)

# View the first few rows
print(df.head())

# Convert values to numeric (just in case)
df = df.apply(pd.to_numeric)

# Get the first row (e.g., B.memory)
first_row = df.iloc[0]

# Sum of the first row
row_sum = first_row.sum()

print(f"Sum of the first row: {row_sum}")