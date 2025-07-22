#usage: python3 nnls.py <normalization>

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
import argparse
from scipy.optimize import nnls
parser = argparse.ArgumentParser()
parser.add_argument("loggedness")
args = parser.parse_args().loggedness

n_samples, n_features = 10, 10321
column_mapping = {
    'B.memory': 'LFQ.intensity.imputed_B.memory_04_steady-state',
    'B.naive': 'LFQ.intensity.imputed_B.naive_04_steady-state',
    'B.plasma': 'LFQ.intensity.imputed_B.plasma_04_steady-state',
    'T4.CM': 'LFQ.intensity.imputed_T4.CM_04_steady-state',
    'T4.EM': 'LFQ.intensity.imputed_T4.EM_04_steady-state',
    'T4.EMRA': 'LFQ.intensity.imputed_T4.EMRA_04_steady-state',
    'T4.naive': 'LFQ.intensity.imputed_T4.naive_04_steady-state',
    'T8.CM': 'LFQ.intensity.imputed_T8.CM_04_steady-state',
    'T8.EM': 'LFQ.intensity.imputed_T8.EM_04_steady-state',
    'T8.EMRA': 'LFQ.intensity.imputed_T8.EMRA_04_steady-state',
    'T8.naive': 'LFQ.intensity.imputed_T8.naive_04_steady-state',
    'Th1': 'LFQ.intensity.imputed_Th1_04_steady-state',
    'Th17': 'LFQ.intensity.imputed_Th17_04_steady-state',
    'Th2': 'LFQ.intensity.imputed_Th2_04_steady-state',
    'mTregs': 'LFQ.intensity.imputed_mTregs_04_steady-state',
    'nTregs': 'LFQ.intensity.imputed_nTregs_04_steady-state',
    'Basophil': 'LFQ.intensity.imputed_Basophil_04_steady-state',
    'Eosinophil': 'LFQ.intensity.imputed_Eosinophil_04_steady-state',
    'MO.classical': 'LFQ.intensity.imputed_MO.classical_04_steady-state',
    'MO.intermediate': 'LFQ.intensity.imputed_MO.intermediate_04_steady-state',
    'MO.nonclassical': 'LFQ.intensity.imputed_MO.nonclassical_04_steady-state',
    'Neutrophil': 'LFQ.intensity.imputed_Neutrophil_04_steady-state',
    'mDC': 'LFQ.intensity.imputed_mDC_04_steady-state',
    'pDC': 'LFQ.intensity.imputed_pDC_04_steady-state',
    'NK.bright': 'LFQ.intensity.imputed_NK.bright_04_steady-state',
    'NK.dim': 'LFQ.intensity.imputed_NK.dim_04_steady-state'
}
fracs_df = pd.read_csv("real_fracs.tsv", sep = "\t", index_col = 0).T
fracs_df_good_cell_order = fracs_df.rename(columns=column_mapping)

sig_matrix = pd.read_csv(f"imputed_sig_matrix_{args}.txt", sep = "\t", index_col=0)
#cell_type_sums = sig_matrix.sum(axis=0)

#sig_matrix = pd.read_csv(f"imputed_sig_matrix_real.txt", sep = "\t", index_col=0) #this was active bevore checking the updated sig_matrix
# sig_matrix = pd.read_csv(f"/Users/samuelgair/Desktop/BA_code/snakemake_fun/BAYESPRISM-updated_sig_matrixdefault_id.txt", sep = "\t", index_col=0).T #this was added
# sig_matrix.reset_index(drop = True, inplace=True)
# print(sig_matrix.shape)
# print(sig_matrix.head)
# sig_matrix.columns = sig_matrix.iloc[0]#until here
# sig_matrix = sig_matrix.multiply(cell_type_sums, axis = 1)
#sig_matrix = sig_matrix[1:]

fracs_df_good_cell_order = fracs_df_good_cell_order[sig_matrix.columns.to_list()] #this was active bevore checking the updated sig_matrix
X = sig_matrix.values


samples_df = pd.read_csv(f"sample_{args}_imputed.txt", sep = "\t", index_col=0)
final_results = []
print("this is samples_df.columns()", samples_df.columns.to_list())
for i in range(len(samples_df.columns.to_list())):
    y = samples_df.iloc[:,i].values#np.dot(X, fracs_np)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5)
    reg_nnls = LinearRegression(positive=True)
    y_pred_nnls = reg_nnls.fit(X_train, y_train).predict(X_test)
    r2_score_nnls = r2_score(y_test, y_pred_nnls)
    coefficients = reg_nnls.coef_ /100
    coefficients = coefficients / coefficients.sum()#added this line
    results = dict(zip(column_mapping.values(), coefficients))
    print("NNLS ", r2_score_nnls)
    print("Coefficients:", coefficients)
    final_results.append(results)

ls_resultsts_df = pd.DataFrame(final_results, index = fracs_df.index.to_list())
ls_resultsts_df.index.name = "Mixture"
ls_resultsts_df.to_csv(f"NNLS-Results_{args}.txt", sep= "\t")
print(f"NNLS-Results_{args}.txt was created")