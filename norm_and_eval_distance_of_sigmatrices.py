import pandas as pd

# Load the matrix (tab-separated, gene IDs as index)
df = pd.read_csv("randomized_bayes_sig_matrix.tsv", sep="\t", index_col=0)

# Convert all values to float (in case some columns are strings)
df = df.apply(pd.to_numeric)

# Normalize each column so it sums to 1
df_normalized = df.div(df.sum(axis=0), axis=1)
print(df_normalized.head())
reference = pd.read_csv("BAYESPRISM-updated_sig_matrix1002.txt", sep="\t", index_col=0)
print(reference.T.head())


normalized_values = normalized.values
reference_values = reference.values

# Compute RMSE per cell type (row-wise)
rmse_per_cell = np.sqrt(np.mean((normalized_values - reference_values) ** 2, axis=1))

# Compute overall RMSE
overall_rmse = np.sqrt(mean_squared_error(reference_values.flatten(), normalized_values.flatten()))

# Print results
print("RMSE per cell type:")
for i, val in enumerate(rmse_per_cell):
    print(f"Cell type {i}: {val:.6f}")

print(f"\nOverall RMSE: {overall_rmse:.6f}")