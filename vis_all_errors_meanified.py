# Example usage: 
# python3 vis_all_errors_meanified.py NNLS-Results_nonlogged.txt NNLS-Results_inlogged.txt NNLS-Results_outlogged.txt CIBERSORT-Results-outlogged.txt CIBERSORT-Results-inlogged.txt CIBERSORT-Results-nonlogged.txt CIBERSORT-Results-QN_nonlogged.txt CIBERSORT-Results-QN_inlogged.txt CIBERSORT-Results-QN_outlogged.txt

import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Define cell groups and column mappings
cell_groups = {
    "B cells": ['B.memory', 'B.naive', 'B.plasma'],
    "T cells": ['T4.CM', 'T4.EM', 'T4.EMRA', 'T4.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'T8.naive', 'Th1', 'Th17', 'Th2', 'mTregs', 'nTregs'],
    "Myeloid cells": ['Basophil', 'Eosinophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'Neutrophil', 'mDC', 'pDC'],
    "NK cells": ['NK.bright', 'NK.dim']
}
column_mapping = {
    'LFQ.intensity.imputed_B.memory_04_steady-state': 'B.memory',
    'LFQ.intensity.imputed_B.naive_04_steady-state': 'B.naive',
    'LFQ.intensity.imputed_B.plasma_04_steady-state': 'B.plasma',
    'LFQ.intensity.imputed_T4.CM_04_steady-state': 'T4.CM',
    'LFQ.intensity.imputed_T4.EM_04_steady-state': 'T4.EM',
    'LFQ.intensity.imputed_T4.EMRA_04_steady-state': 'T4.EMRA',
    'LFQ.intensity.imputed_T4.naive_04_steady-state': 'T4.naive',
    'LFQ.intensity.imputed_T8.CM_04_steady-state': 'T8.CM',
    'LFQ.intensity.imputed_T8.EM_04_steady-state': 'T8.EM',
    'LFQ.intensity.imputed_T8.EMRA_04_steady-state': 'T8.EMRA',
    'LFQ.intensity.imputed_T8.naive_04_steady-state': 'T8.naive',
    'LFQ.intensity.imputed_Th1_04_steady-state': 'Th1',
    'LFQ.intensity.imputed_Th17_04_steady-state': 'Th17',
    'LFQ.intensity.imputed_Th2_04_steady-state': 'Th2',
    'LFQ.intensity.imputed_mTregs_04_steady-state': 'mTregs',
    'LFQ.intensity.imputed_nTregs_04_steady-state': 'nTregs',
    'LFQ.intensity.imputed_Basophil_04_steady-state': 'Basophil',
    'LFQ.intensity.imputed_Eosinophil_04_steady-state': 'Eosinophil',
    'LFQ.intensity.imputed_MO.classical_04_steady-state': 'MO.classical',
    'LFQ.intensity.imputed_MO.intermediate_04_steady-state': 'MO.intermediate',
    'LFQ.intensity.imputed_MO.nonclassical_04_steady-state': 'MO.nonclassical',
    'LFQ.intensity.imputed_Neutrophil_04_steady-state': 'Neutrophil',
    'LFQ.intensity.imputed_mDC_04_steady-state': 'mDC',
    'LFQ.intensity.imputed_pDC_04_steady-state': 'pDC',
    'LFQ.intensity.imputed_NK.bright_04_steady-state': 'NK.bright',
    'LFQ.intensity.imputed_NK.dim_04_steady-state': 'NK.dim'
}

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("errors", nargs="*")  # Input error files
args = parser.parse_args()

# Load results
results = {}
for resultfile in args.errors:
    result_name = os.path.splitext(os.path.basename(resultfile))[0]
    results[result_name] = pd.read_csv(resultfile, sep="\t")

# Load real fractions
real_fracs_coarse_compare = pd.read_csv("real_coarse_fracs.tsv", sep='\t').drop(columns=["Unnamed: 0"])
real_fracs_coarse_compare.index = ["NK cells, real", "B cells, real", "Myeloid cells, real", "T cells, real"]

real_fracs_compare = pd.read_csv("real_fracs.tsv", sep='\t')
real_fracs_compare.index = real_fracs_compare["Unnamed: 0"]
real_fracs_compare = real_fracs_compare.drop(columns=["Unnamed: 0"]).T
real_fracs_compare.index = map(lambda a: a + "_predict", real_fracs_compare.index.tolist())

# Calculate average absolute error for each method
method_avg_errors = {}
for key in results.keys():
    results_cibersort = results[key].rename(columns=column_mapping)
    combined_coarse_df = pd.concat([results_cibersort.reset_index(drop=True), real_fracs_coarse_compare.T.reset_index(drop=True)], axis=1)

    # Sum the predicted subtypes into coarse types
    for coarse_type, subtypes in cell_groups.items():
        combined_coarse_df[coarse_type] = combined_coarse_df[subtypes].sum(axis=1)
    combined_coarse_df = combined_coarse_df.loc[:, combined_coarse_df.columns.str.contains("cells")]
    predicted_to_real_name = {key: f"{key}, real" for key in cell_groups.keys()}

    error_data = {cell_type: combined_coarse_df[cell_type] - combined_coarse_df[predicted_to_real_name[cell_type]] for cell_type in cell_groups.keys()}
    error_df = pd.DataFrame(error_data)

    # Calculate the average absolute error for the current method
    values = error_df.values.flatten()
    method_avg_error = np.mean(np.abs(values))
    method_avg_errors[key] = method_avg_error

# Sort methods by average absolute error
sorted_methods = sorted(method_avg_errors.keys(), key=lambda x: method_avg_errors[x])

# Plot the errors
plt.figure(figsize=(10, 6))
color_mapping = {
    "B cells": "red",
    "T cells": "blue",
    "Myeloid cells": "orange",
    "NK cells": "black"
}
cell_type_offsets = {
    "B cells": -0.15,
    "T cells": -0.05,
    "Myeloid cells": 0.05,
    "NK cells": 0.15
}
#plt.tick_params(axis='y', labelsize = 5)

for i, key in enumerate(sorted_methods):
    results_cibersort = results[key].rename(columns=column_mapping)
    combined_coarse_df = pd.concat([results_cibersort.reset_index(drop=True), real_fracs_coarse_compare.T.reset_index(drop=True)], axis=1)

    # Sum the predicted subtypes into coarse types
    for coarse_type, subtypes in cell_groups.items():
        combined_coarse_df[coarse_type] = combined_coarse_df[subtypes].sum(axis=1)
    combined_coarse_df = combined_coarse_df.loc[:, combined_coarse_df.columns.str.contains("cells")]
    predicted_to_real_name = {key: f"{key}, real" for key in cell_groups.keys()}

    error_data = {cell_type: combined_coarse_df[cell_type] - combined_coarse_df[predicted_to_real_name[cell_type]] for cell_type in cell_groups.keys()}
    error_df = pd.DataFrame(error_data)

    values = error_df.values.flatten()
    cell_types = np.repeat(error_df.columns, error_df.shape[0])
    colors = [color_mapping[cell] for cell in cell_types]

    # Plot scatter points for this method
    x_values = np.array([i + cell_type_offsets[cell_type] for cell_type in cell_types])
    plt.scatter(x=x_values, y=values, c=colors, edgecolor="black", s=10, alpha=0.8)

    # Bar plot for average error per cell type
    cell_type_avg_errors = error_df.mean(axis=0)
    for cell_type, avg_error in cell_type_avg_errors.items():
        x_value = i + cell_type_offsets[cell_type]
        plt.bar(
            x=x_value,
            height=avg_error,
            width=0.1,
            color=color_mapping[cell_type],
            alpha=0.7,
            label=f"{cell_type} Avg Error" if i == 0 else None
        )

    # Add a horizontal line for the average absolute error of the current method
    method_avg_error = np.mean(np.abs(values))
    plt.axhline(
        y=method_avg_error,
        color="green",
        linestyle="--",
        xmin=i / len(sorted_methods),
        xmax=(i + 1) / len(sorted_methods)
    )
    plt.text(
        x=i + 0.15,
        y=method_avg_error + 0.01,
        s=f"{method_avg_error:.4f}",
        color="green",
        fontsize=15,
        ha="left",
        va="center"
    )

# Add legend
for cell_type, color in color_mapping.items():
    plt.scatter([], [], color=color, label=cell_type)
plt.scatter([], [], color="green", marker="_", label="Average Absolute Errors")
plt.axhline(y=0, color="black", linestyle="-", linewidth=1)

plt.xticks(range(len(sorted_methods)), sorted_methods, fontsize=12, rotation=45)
plt.xlabel("Deconvolution Methods (Sorted by Avg Abs Error)", fontsize=10)
plt.ylabel("Difference between predicted fractions and real fractions")
plt.title("Error Scatter Plot by Cell Type and Method")
plt.legend(loc="lower right", fontsize = 12)
plt.tight_layout()

# Show the plot
plt.show()