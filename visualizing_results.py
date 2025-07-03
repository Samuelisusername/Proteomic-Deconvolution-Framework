#example usage can be used for when 
import pandas as pd
import numpy as np
import argparse
parser = argparse.ArgumentParser(description= "Vizualize prediction results.")
parser.add_argument("results")
parser.add_argument("method")


args = parser.parse_args()
results = pd.read_csv(args.results, sep = '\t', index_col = 0)


#cibersort results
# results_cibersort = pd.read_csv("CIBERSORT-Results.txt", sep = '\t')
# results_cibersort.index = results_cibersort["Mixture"]
# results_cibersort = results_cibersort.drop(columns=["Mixture"])
# results_cibersort = results_cibersort.drop(columns=["Correlation"])
# results_cibersort = results_cibersort.drop(columns=["RMSE"])
# results_cibersort = results_cibersort.drop(columns=["P-value"]) 

results_cibersort = results #changed for snakemake

#coarse real fracs
real_fracs_coarse_compare = pd.read_csv("real_coarse_fracs.tsv", sep = '\t') #should it not be txt??
real_fracs_coarse_compare = real_fracs_coarse_compare.drop(columns = ["Unnamed: 0"])
real_fracs_coarse_compare.index = ["NK cells, real", "B cells, real", "Myeloid cells, real","T cells, real"] #i really hope this is correct..


#fine real fracs
real_fracs_compare = pd.read_csv("real_fracs.tsv", sep = '\t')
real_fracs_compare.index = real_fracs_compare["Unnamed: 0"]
real_fracs_compare= real_fracs_compare.drop(columns = ["Unnamed: 0"]).T
real_fracs_compare.index = map(lambda a: a + "_predict", real_fracs_compare.index.tolist())

# Rename the columns in results_cibersort to match real_fracs_compare
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
results_cibersort = results_cibersort.rename(columns=column_mapping) #this has no effect when doing the coarse

#fine df
combined_df = pd.concat([results_cibersort, real_fracs_compare], axis=0)
print(combined_df)
#combined_coarse_df = pd.concat([results_cibersort, real_fracs_coarse_compare.T], axis = 0)

#coarse df
combined_coarse_df = pd.concat([results_cibersort.reset_index(drop=True), real_fracs_coarse_compare.T.reset_index(drop=True)], axis=1)



cell_groups = {
    "B cells": ['B.memory', 'B.naive', 'B.plasma'],
    "T cells": ['T4.CM', 'T4.EM', 'T4.EMRA', 'T4.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'T8.naive', 'Th1', 'Th17', 'Th2', 'mTregs', 'nTregs'],
    "Myeloid cells": ['Basophil', 'Eosinophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'Neutrophil', 'mDC', 'pDC'],
    "NK cells": ['NK.bright', 'NK.dim']
}

# Calculate sums for each group
group_sums = {}
for group_name, columns in cell_groups.items():
    group_sums[group_name] = combined_df[columns].sum(axis=1) #summing up the fracs to get the coarse type. this is real sums, not predicred

# Create a new DataFrame for the group sums
group_sums_df = pd.DataFrame(group_sums)


import matplotlib.pyplot as plt

# Transpose the DataFrame for easier plotting
combined_df_T = combined_df.T
group_sums_df_T = group_sums_df.T

""" # Plot for each sample
combined_coarse_df.T.plot(kind='bar', figsize=(15, 7), width=0.8)
plt.title("real and predicted ones")
plt.ylabel("Proportion")
plt.xlabel("Cell Types")
plt.legend(loc="upper right")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()  # Save as PNG """
#print(combined_coarse_df) """

""" group_sums_df_T.plot(kind='bar', figsize=(12, 6))
plt.title("Sum of Cell Type Groups for Each Sample")
plt.ylabel("Proportion")
plt.xlabel("Samples")
plt.legend(title="Cell Type Groups")
plt.legend().remove()
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("cell_type_group_sums.png")  # Save as PNG """

def visualize_predicted_vs_real(df):
    global method
    """
    Generates stacked bar charts comparing predicted and real cell fractions.

    Parameters:
        df (pd.DataFrame): DataFrame with the following columns:
            - 'B cells', 'T cells', 'Myeloid cells', 'NK cells'
            - 'B Cells real', 'T cells real', 'Myeloid cells real', 'NK cells real'
    """
    # Split into predicted and real
    #predicted_cols = ['B cells', 'T cells', 'Myeloid cells', 'NK cells']
    cell_groups = {
        "B cells": ['B.memory', 'B.naive', 'B.plasma'],
        "T cells": ['T4.CM', 'T4.EM', 'T4.EMRA', 'T4.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'T8.naive', 'Th1', 'Th17', 'Th2', 'mTregs', 'nTregs'],
        "Myeloid cells": ['Basophil', 'Eosinophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'Neutrophil', 'mDC', 'pDC'],
        "NK cells": ['NK.bright', 'NK.dim']
    }

    # Sum the predicted subtypes into coarse types
    for coarse_type, subtypes in cell_groups.items():
        df[coarse_type] = df[subtypes].sum(axis=1)
    df = df.loc[:, df.columns.str.contains("cells")]
    print(df)

    # Define predicted and real coarse type columns
    predicted_cols = ['B cells', 'T cells', 'Myeloid cells', 'NK cells']
    real_cols = ['B cells, real', 'T cells, real', 'Myeloid cells, real', 'NK cells, real']
    #print(df.columns)
    fig, axs = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    # Plot predicted values
    df[predicted_cols].plot(kind='bar', stacked=True, ax=axs[0], colormap='tab20')
    axs[0].set_title("Predicted Cell Fractions")
    axs[0].set_ylabel("Fraction")
    axs[0].legend(loc="upper right")

    # Plot real values
    df[predicted_cols+real_cols].to_csv(f"{args.method}_fractions_compare.csv")
    df[real_cols].plot(kind='bar', stacked=True, ax=axs[1], colormap='tab20')
    axs[1].set_title("Real Cell Fractions")
    axs[1].set_xlabel("Sample")
    axs[1].set_ylabel("Fraction")
    axs[1].legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(f"all_fracs_compare/{args.method}_fractions_compare.png")
def visualize_prediciton_errors(df):
    global real_fracs_coarse_compare
    cell_groups = {
    "B cells": ['B.memory', 'B.naive', 'B.plasma'],
    "T cells": ['T4.CM', 'T4.EM', 'T4.EMRA', 'T4.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'T8.naive', 'Th1', 'Th17', 'Th2', 'mTregs', 'nTregs'],
    "Myeloid cells": ['Basophil', 'Eosinophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'Neutrophil', 'mDC', 'pDC'],
    "NK cells": ['NK.bright', 'NK.dim']
    }

    # Sum the predicted subtypes into coarse types
    for coarse_type, subtypes in cell_groups.items():
        df[coarse_type] = df[subtypes].sum(axis=1)
    df = df.loc[:, df.columns.str.contains("cells")]
    predicted_to_real_name = {key: f"{key}, real" for key in cell_groups.keys()}

    error_data = {cell_type:df[cell_type]-df[predicted_to_real_name[cell_type]] for cell_type in cell_groups.keys()}
    avg_error = [(abs(df[cell_type]-df[predicted_to_real_name[cell_type]])).mean() for cell_type in cell_groups.keys()]
    print("this is avg_error : \n",np.mean(avg_error))
    error_df = pd.DataFrame(error_data)


    # Plot
    error_df.plot(kind='bar', figsize=(14, 6))
    error_df.to_csv("Errors_of_coarse_celltype_prediction_Cibersort.csv")
    plt.title("Error Between Predicted and Real Cell Fractions")
    plt.xlabel("Sample")
    plt.ylabel("Error")
    plt.legend(title="Cell Type", loc="upper right")
    plt.tight_layout()
    plt.savefig(f"all_errors/{args.method}_errors.png")
    plt.close()
    for celltype in cell_groups.keys():
        plt.scatter(df[predicted_to_real_name[celltype]], abs(error_df[celltype]), label = celltype)
    plt.xlabel('real percentage')  # Replace with your x-axis label
    plt.ylabel('error')  # Replace with your y-axis label
    plt.title('Scatter Plot of percent vs error')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"correlation_images/{args.method}_all_cells_percents_error_correlation.png")

def visualize_predicted_vs_real_percentages(df):
    cell_groups = {
    "B cells": ['B.memory', 'B.naive', 'B.plasma'],
    "T cells": ['T4.CM', 'T4.EM', 'T4.EMRA', 'T4.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'T8.naive', 'Th1', 'Th17', 'Th2', 'mTregs', 'nTregs'],
    "Myeloid cells": ['Basophil', 'Eosinophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'Neutrophil', 'mDC', 'pDC'],
    "NK cells": ['NK.bright', 'NK.dim']}
    # Sum the predicted subtypes into coarse types
    for coarse_type, subtypes in cell_groups.items():
        df[coarse_type] = df[subtypes].sum(axis=1)
    df = df.loc[:, df.columns.str.contains("cells")]
    predicted_to_real_name = {key: f"{key}, real" for key in cell_groups.keys()}

    for celltype in cell_groups.keys():
        new_df = df[[col for col in df.columns if celltype in col]]
        plt.figure(figsize=(8, 6))
        plt.scatter(new_df[celltype], new_df[predicted_to_real_name[celltype]], label=celltype)
        plt.xlabel('predicted')  # Replace with your x-axis label
        plt.ylabel('real')  # Replace with your y-axis label
        plt.title('Scatter Plot of X vs Y')
        plt.legend()
        plt.grid(True)
        plt.savefig(f"correlation_images/{args.method}_{celltype}_percents_correlation.png")

        plt.figure(figsize=(10, 8))  # Create a single figure for all plots

    for celltype in cell_groups.keys():
        new_df = df[[col for col in df.columns if celltype in col]]
        plt.scatter(
            new_df[celltype], 
            new_df[predicted_to_real_name[celltype]], 
            label=celltype  # Use celltype as the label for the legend
        )

    # Add labels, title, legend, and grid
    plt.xlabel('Predicted')  # Replace with your x-axis label
    plt.ylabel('Real')  # Replace with your y-axis label
    plt.title('Scatter Plot of Predicted vs Real for All Cell Types')
    plt.legend(title="Cell Types")  # Add a legend with a title
    plt.grid(True)

    # Save the combined plot
    plt.savefig(f"correlation_images/{args.method}_all_celltypes_percents_correlation.png")
# def visualize_predicted_percentage_vs_error():





#visualize_predicted_vs_real(combined_coarse_df)
#print(combined_coarse_df)
#visualize_prediciton_errors(combined_coarse_df)
""" errors = pd.DataFrame({
    'B cells': combined_coarse_df['B cells'] - combined_coarse_df['B cells, real'],
    'T cells': combined_coarse_df['T cells'] - combined_coarse_df['T cells, real'],
    'Myeloid cells': combined_coarse_df['Myeloid cells'] - combined_coarse_df['Myeloid cells, real'],
    'NK cells': combined_coarse_df['NK cells'] - combined_coarse_df['NK cells, real'],
})
errorsdic = {
    'B cells': combined_coarse_df['B cells'] - combined_coarse_df['B cells, real'],
    'T cells': combined_coarse_df['T cells'] - combined_coarse_df['T cells, real'],
    'Myeloid cells': combined_coarse_df['Myeloid cells'] - combined_coarse_df['Myeloid cells, real'],
    'NK cells': combined_coarse_df['NK cells'] - combined_coarse_df['NK cells, real'],
}
abserrorlist = [abs(errorsdic[something])for something in errorsdic.keys()]
print("mean error of coarse data: \n", np.mean(abserrorlist))
errors.plot(kind = 'bar')
plt.title("Prediction Errors")
plt.show() """
visualize_prediciton_errors(combined_coarse_df)
visualize_predicted_vs_real(combined_coarse_df)
visualize_predicted_vs_real_percentages(combined_coarse_df)
