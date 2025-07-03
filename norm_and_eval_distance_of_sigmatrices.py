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