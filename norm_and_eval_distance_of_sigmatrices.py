# usage python3 norm_and_eval_distance_of_sigmatrices <real_signature_matrix> <reference_signature_matrix> <real_signature_matrix>
# example python3 norm_and_eval_distance_of_sigmatrices.py "real_sig_nonlogged.txt" "randomized_bayes_sig_matrix.tsv" "BAYESPRISM-updated_sig_matrixdefault_id.txt"

import matplotlib.pyplot
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
import os 
import argparse
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Evaluate distance of signature matrices")
parser.add_argument("real_sig", help="Path to the real signature matrix")
parser.add_argument("ref_sig", help="Path to the reference signature matrix")
parser.add_argument("updated_sig", help="Path to the updated signature matrix")
args = parser.parse_args()

# Load the matrix (tab-separated, gene IDs as index)
df = pd.read_csv(args.ref_sig, sep="\t", index_col=0)
df_normalized_ref = df.div(df.sum(axis=0), axis=1).T
ref_cols = list(df_normalized_ref.T.columns)


real = pd.read_csv(args.real_sig, index_col = 0, sep = "\t")
df_normalized_real = real.div(real.sum(axis = 0), axis = 1).T

#print(df_normalized_ref.head())
updated_reference = pd.read_csv(args.updated_sig, sep="\t", index_col=0).T
updated_cols = list(updated_reference.columns)
print(updated_cols, "this was updated_cols")
print(len(list(updated_reference.index)), " this is len of index")

df_normalized_real = df_normalized_real.T
df_normalized_ref = df_normalized_ref.T

print(df_normalized_real.head)
print(df_normalized_ref.head, "first")
desired_order = [
    'B.memory', 'B.naive', 'B.plasma', 'mTregs', 'T4.naive', 'nTregs', 'T4.CM', 'T4.EM', 'T4.EMRA',
    'Th1', 'Th17', 'Th2', 'T8.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'mDC', 'pDC', 'Basophil',
    'Eosinophil', 'Neutrophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'NK.bright', 'NK.dim'
]
df_normalized_real = df_normalized_real[desired_order]
print(df_normalized_real.head, "second")
df_normalized_real.reset_index(drop=True,inplace = True)
df_normalized_ref.reset_index(drop = True, inplace = True)
updated_reference.reset_index(drop= True, inplace= True)

df_normalized_real.columns = [f"{col}-Real" for col in df_normalized_real.columns]
df_normalized_ref.columns = [f"{col}-Ref" for col in df_normalized_ref.columns]
updated_reference.columns = [f"{col}-updated" for col in updated_reference.columns]

combined_df = pd.concat([df_normalized_real, df_normalized_ref, updated_reference], axis="columns")
print(len(combined_df.index), "this is the combined inndex length")
print(combined_df.head)
#print(combined_df)
pca = PCA(n_components=2)
pca_results = pca.fit_transform(combined_df.T)
print(pca_results[:5])
pca_df = pd.DataFrame(pca_results, columns = ["PC1", "PC2"])
print(pca_df.shape)
print(pca_df.index, "this is pca index")
pca_df["source"]= ["Real"] * len(df_normalized_real.columns) + ["Ref"] * len(df_normalized_ref.columns) + ["Updated"] * len(updated_reference.columns)
#print(len(pca_df["source"]))
print(len(list(real.columns)), len(ref_cols), len(updated_cols))
pca_df["Cell type"] =  updated_cols * 3

pd.set_option("display.max_rows", None)
print(pca_df)
plt.figure()
colours = {"Real": "blue", "Ref": "green", "Updated" : "red"}
for source in colours.keys():
    subset = pca_df[pca_df["source"] == source]
    plt.scatter(subset["PC1"], subset["PC2"], label = source, color = colours[source], alpha = 0.7) 

    
ref_points = pca_df[pca_df["source"] == "Ref"]
updated_points = pca_df[pca_df["source"] == "Updated"]
real_points = pca_df[pca_df["source"]=="Real"]

for i in range(len(ref_points)):
    # Green arrows from Ref to Real
    plt.quiver(
        ref_points.iloc[i]["PC1"], ref_points.iloc[i]["PC2"],  # Start point (Ref)
        real_points.iloc[i]["PC1"] - ref_points.iloc[i]["PC1"],  # Delta X (Real - Ref)
        real_points.iloc[i]["PC2"] - ref_points.iloc[i]["PC2"],  # Delta Y (Real - Ref)
        angles='xy', scale_units='xy', scale=1, color="green", alpha=0.7, width = 0.001
    )
    # Red arrows from Ref to Updated
    plt.quiver(
        ref_points.iloc[i]["PC1"], ref_points.iloc[i]["PC2"],  # Start point (Ref)
        updated_points.iloc[i]["PC1"] - ref_points.iloc[i]["PC1"],  # Delta X (Updated - Ref)
        updated_points.iloc[i]["PC2"] - ref_points.iloc[i]["PC2"],  # Delta Y (Updated - Ref)
        angles='xy', scale_units='xy', scale=1, color="red", alpha=0.7, width = 0.001
    )


# for i, row in pca_df.iterrows():
#     plt.annotate(row["Cell type"], (row["PC1"], row["PC2"]), fontsize=6, alpha=0.7)
from adjustText import adjust_text

texts = []
for i, row in pca_df[pca_df["source"] == "Real"].iterrows():
    texts.append(
        plt.text(row["PC1"], row["PC2"], row["Cell type"], fontsize=10, alpha=0.8)
    )

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.title("PCA Scatter Plot of Signature Matrices")
plt.legend(fontsize = 15)
plt.grid()

# Show the plot
plt.show()






normalized_values_ref = df_normalized_ref.values
normalized_real_values= df_normalized_real.values

# Compute RMSE per cell type (row-wise)
# rmse_per_cell = np.sqrt(np.mean((normalized_values_ref - reference_values) ** 2, axis=1))

# Compute overall RMSE
#overall_rmse = np.sqrt(mean_squared_error(normalized_real_values.flatten(), normalized_values_ref.flatten()))
real_v_reference = np.sqrt(mean_squared_error(normalized_real_values.flatten(), normalized_values_ref.flatten()))
real_v_updated = np.sqrt(mean_squared_error(normalized_real_values.flatten(), updated_reference.values.flatten()))
print("real v ref = ", real_v_reference, "; real v updated =   ", real_v_updated)
