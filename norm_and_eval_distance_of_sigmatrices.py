# usage python3 norm_and_eval_distance_of_sigmatrices "real_sig_nonlogged.txt" "randomized_bayes_sig_matrix.tsv" "BAYESPRISM-updated_sig_matrixdefault_id.txt"

import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
import os 
import argparse

parser = argparse.ArgumentParser(description="Evaluate distance of signature matrices")
parser.add_argument("real_sig", help="Path to the real signature matrix")
parser.add_argument("ref_sig", help="Path to the reference signature matrix")
parser.add_argument("updated_sig", help="Path to the updated signature matrix")
args = parser.parse_args()

# Load the matrix (tab-separated, gene IDs as index)
df = pd.read_csv(args.ref_sig, sep="\t", index_col=0)
df_normalized_ref = df.div(df.sum(axis=0), axis=1)

real = pd.read_csv(args.real_sig, index_col = 0, sep = "\t")
df_normalized_real = real.div(real.sum(axis = 0), axis = 1)

print(df_normalized_ref.head())
updated_reference = pd.read_csv(args.updated_sig, sep="\t", index_col=0) #you gotta manually adapt this name

df_normalized_real = df_normalized_real.T
df_normalized_ref = df_normalized_ref.T


normalized_values_ref = df_normalized_ref.values
normalized_real_values= df_normalized_real.values

# Compute RMSE per cell type (row-wise)
# rmse_per_cell = np.sqrt(np.mean((normalized_values_ref - reference_values) ** 2, axis=1))

# Compute overall RMSE
#overall_rmse = np.sqrt(mean_squared_error(normalized_real_values.flatten(), normalized_values_ref.flatten()))
real_v_reference = np.sqrt(mean_squared_error(normalized_real_values.flatten(), normalized_values_ref.flatten()))
real_v_updated = np.sqrt(mean_squared_error(normalized_real_values.flatten(), updated_reference.values.flatten()))
print("real v ref = ", real_v_reference, "; real v updated =   ", real_v_updated)