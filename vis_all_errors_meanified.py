#example usage: 
#python3 vis_all_errors_meanified.py NNLS-Results_nonlogged.txt
import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

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
parser = argparse.ArgumentParser()
parser.add_argument("errors", nargs="*")#the idea of this file is that for the same sample and same sig_matrix we compare all the methods
args = parser.parse_args()
results = {}
for i, resultfile in enumerate(args.errors):
    result_name = os.path.splitext(os.path.basename(resultfile))[0] 
    results[result_name] = pd.read_csv(resultfile, sep="\t")


real_fracs_coarse_compare = pd.read_csv("real_coarse_fracs.tsv", sep = '\t') #should it not be txt??
real_fracs_coarse_compare = real_fracs_coarse_compare.drop(columns = ["Unnamed: 0"])
real_fracs_coarse_compare.index = ["NK cells, real", "B cells, real", "Myeloid cells, real","T cells, real"] #i really hope this is correct..


#fine real fracs
real_fracs_compare = pd.read_csv("real_fracs.tsv", sep = '\t')
real_fracs_compare.index = real_fracs_compare["Unnamed: 0"]
real_fracs_compare= real_fracs_compare.drop(columns = ["Unnamed: 0"]).T
real_fracs_compare.index = map(lambda a: a + "_predict", real_fracs_compare.index.tolist())


for i,key in enumerate(results):
    print(i, "this is i")
    cell_groups = {
    "B cells": ['B.memory', 'B.naive', 'B.plasma'],
    "T cells": ['T4.CM', 'T4.EM', 'T4.EMRA', 'T4.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'T8.naive', 'Th1', 'Th17', 'Th2', 'mTregs', 'nTregs'],
    "Myeloid cells": ['Basophil', 'Eosinophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'Neutrophil', 'mDC', 'pDC'],
    "NK cells": ['NK.bright', 'NK.dim']
    }
    results_cibersort = results[key]
    results_cibersort = results_cibersort.rename(columns=column_mapping) #this has no effect when doing the coarse
    combined_df = pd.concat([results_cibersort, real_fracs_compare], axis=0)
    combined_coarse_df = pd.concat([results_cibersort.reset_index(drop=True), real_fracs_coarse_compare.T.reset_index(drop=True)], axis=1)
    group_sums = {}
    for group_name, columns in cell_groups.items():
        group_sums[group_name] = combined_df[columns].sum(axis=1) #summing up the fracs to get the coarse type. this is real sums, not predicred
    group_sums_df = pd.DataFrame(group_sums)

    # Sum the predicted subtypes into coarse types
    for coarse_type, subtypes in cell_groups.items():
        combined_coarse_df[coarse_type] = combined_coarse_df[subtypes].sum(axis=1)
    combined_coarse_df = combined_coarse_df.loc[:, combined_coarse_df.columns.str.contains("cells")]
    predicted_to_real_name = {key: f"{key}, real" for key in cell_groups.keys()}

    error_data = {cell_type:combined_coarse_df[cell_type]-combined_coarse_df[predicted_to_real_name[cell_type]] for cell_type in cell_groups.keys()}
    avg_error = [(abs(combined_coarse_df[cell_type]-combined_coarse_df[predicted_to_real_name[cell_type]])).mean() for cell_type in cell_groups.keys()]
    print()
    print("this is avg_error : \n",np.mean(avg_error))
    error_df = pd.DataFrame(error_data)
    print(error_df)

    values = error_df.values.flatten()
    print(values, "this was values")
    cell_types = np.repeat(error_df.columns, error_df.shape[0])  # Repeat column names for each value
    #new way of repreating

    cell_types = np.concatenate([error_df.columns, error_df.columns, error_df.columns, error_df.columns, error_df.columns, error_df.columns, error_df.columns, error_df.columns, error_df.columns, error_df.columns])
    print(cell_types, "this is cell_types")

    # Assign colors to each cell type
    color_mapping = {
        "B cells": "red",
        "T cells": "blue",
        "Myeloid cells": "orange",
        "NK cells": "black"
    }
    colors = [color_mapping[cell] for cell in cell_types]


    # Calculate the overall average error
    overall_avg_error = np.mean(values)

    cell_type_offsets = {
    "B cells": -0.15,  # Offset for B cells
    "T cells": -0.05,  # Offset for T cells
    "Myeloid cells": 0.05,  # Offset for Myeloid cells
    "NK cells": 0.15  # Offset for NK cells
    }
    # Plot the scatter plot
    x_values = np.array([i + cell_type_offsets[cell_type] for cell_type in cell_types])  # Set x-values to the current iteration index

    # Plot the scatter points for this iteration

    plt.scatter(x=x_values, y=values, c=colors, edgecolor="black", s=10, alpha=0.8)

    cell_type_avg_errors = error_df.mean(axis=0)  # Mean error for each cell type
    cell_type_x_values = [i] * len(cell_type_avg_errors) 

    for cell_type, avg_error in cell_type_avg_errors.items():
        group_vals = error_df[cell_type]
        x_value = i + cell_type_offsets[cell_type]
        y_value = group_vals.mean()
        label = f"{cell_type}Avg Error" if i==0 else None
        plt.bar(
            x = x_value,
            #x=cell_type_x_values[cell_type_avg_errors.index.get_loc(cell_type)],
            #y=y_value,
            height= y_value,
            color=color_mapping[cell_type],
            #marker="*",  # Star marker
            width= 0.1,
            alpha = 0.7,
            #s=100,  # Larger size for visibility
            label=label
        )
        print(f"{key} - B cells min/max: {error_df['B cells'].min():.4f} / {error_df['B cells'].max():.4f}")
        print(f"{key} - B cells mean: {error_df['B cells'].mean():.4f}")

        # Calculate the average error for the current method
        method_avg_error = np.mean(values)
        print(f"{key} - Method Avg Error: {method_avg_error:.4f}")

        # Add a horizontal line for the average error of the current method
        # plt.axhline(
        #     y=method_avg_error,
        #     color="green",
        #     linestyle="--",
        #     #label=f"{key} Avg Error: {method_avg_error:.4f}",
        #     xmin=i / len(results.keys()),  # Start of the line for this method
        #     xmax=(i + 1) / len(results.keys())  # End of the line for this method
        # )


        # Calculate the absolute average error for the current method
        print(values, "this is values")
        method_avg_error = np.mean(np.abs(values))  # Take the absolute value of errors before calculating the mean


        # Add a horizontal line for the absolute average error of the current method
        plt.axhline(
            y=method_avg_error,
            color="green",
            linestyle="--",
            #label=f"{key} Abs Avg Error: {method_avg_error:.4f}",
            xmin=i / len(results.keys()),  # Start of the line for this method
            xmax=(i + 1) / len(results.keys())  # End of the line for this method
        )

        plt.text(
            x=i + 0.1,  # Position the text near the middle of the method's range
            y=method_avg_error + 0.01,  # Slightly above the line
            s=f"{method_avg_error:.4f}",  # Display the average absolute error value
            color="green",
            fontsize=8,
            ha="left",  # Center-align the text
            va = "center"
        )
    if i==0:
        plt.scatter([], [], color="green", marker="_", label="Average Absolute Errors")

# Add a horizontal line for the overall average error
# plt.axhline(y=overall_avg_error, color="green", linestyle="--")

for cell_type, color in color_mapping.items():
    plt.scatter([], [], color=color, label=cell_type)  # Empty scatter for legend
plt.legend(loc="upper right", title="Deconvolution Methods")
# Add labels and legend
plt.xticks(range(len(results.keys())), results.keys(), fontsize=6, rotation = 45)
plt.xlabel("Deconvolution Methods", fontsize = 6)
plt.ylabel("Difference between predicted fractions and real fractions")
plt.title("Error Scatter Plot by Cell Type and Method")
plt.legend(loc="lower right")
plt.scatter([], [], color="green", marker="_", label="Average Absolute Errors for each method")

plt.axhline(
    y=0,
    color="black",
    linestyle="-",  # Solid line
    linewidth=1,    # Line thickness
    label="y = 0"   # Add a label for the legend
)
plt.tight_layout()


# Show the plot
plt.show()