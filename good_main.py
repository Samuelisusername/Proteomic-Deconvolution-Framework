#usage:
    # python3 good_main.py <number_of_healthy_bulk_samples> <number_of_cancerous_bulk_samples> <method> <job_id>

#example usage:
    # python3 good_main.py 5 5 ciber_and_nnls 50

#overview: 
# This script generates data for deconvolution method benchmarking purposes. 
# To that end it provides the user with a signature matrix as an input to the method to be tested,  samples that the method can be tested on , and the ground truth cell composition of the samples. 
# (if the bayes_or_ciber_and_nnls parameter is set as bayes it returns myinput2.gbm.rdata ontop of the previously mentioned outputs. This rdata is formatted such that Bayesprism can be run on those samples)

#1. <normalization>
        # The user decides what kind of normalization they want to use for the data: 
            # inlogged: *Bevore* any proteomic profiles are combined to a bulk sample using weighted averages they undergo log transformation (proteomic_profile -> np.log1p(proteomic_profile)
            # outlogged: *After* all proteomic profiles are combined to a bulk sample using weighted averages they undergo log transformation (proteomic_profile -> np.log1p(proteomic_profile)
            # nonlogged: The Protein intensities are never logged.
#for ciber and nnls we have all three. 
#For bayesprism atm only nonlogged


# input parameters: 

    #1&2. <number_of_healthy_bulk_samples> <number_of_cancerous_bulk_samples>
        # The user decides on the number of healthy and cancerous samples that they want to generate
    #3. <method>
        # The user decides if the output data should be made for cibersort and other methods that require .txt files, or they can decide on curating the output the Bayesprism method (write "ciber_and_nnls" or "bayes" for your preferred method)
        #(Note: if you dont know which to use or just want samples to adjust for whatever fancy other methods you want to use, just use "ciber_and_nnls", "bayes" is really specific to bayesprism and you'll probs have a hard time altering that format)
    #4. <job_id>
        #this is to help the user assign each iteration of running this code a unique number.  for example: in case you want to test the method on different signature matrices, then you just run this code several times with different job_id's and your signature matrices will be saved in the sigs directory with that specific job id.
        #Usually you wont need this, and in that case just write some jibberish integer to fill the spot
# outputs: 
    #imputed_sig_matrix_<normalization>.txt
        #Cibersort and NNLS(Non-negative least squares) read from this and treat it as they proteomic signature matrix
    #healthy_sample_imputed_<normalization>.txt
        #Cibersort and NNLS use this as their sample file. don't be confused about the "healthy" word in the beginning of the file name: It will contain healthy and cancrous samples, as many as you told the code to generate (<number_of_healthy_bulk_samples> <number_of_cancerous_bulk_samples>). The healthy stands for the fact that the reference proein profiles are from healthy patients. 
    #real_fracs.tsv, real_coarse_fracs.tsv
        #This are the ground truth fractions of the samples you've given. once the ground truth fraction for all 26 celltypes (real_fracs.tsv) and once for the fractions of the 4 main celltypes that the 26 celltypes aggregate to. (For example: B.memory + B.naive + B.plasma = B Cells)
    # there are some more outputs: those are for debugging purposes or for the purpose of saving sigmatrices across runs. 

#How it works and other questions: 
    #read my thesis or contact me at sagair@ethz.ch




import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import sys


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("nr_of_healthy_samples")
parser.add_argument("nr_unhealthy")
#parser.add_argument("randsig?")
parser.add_argument("bayes_or_ciber_and_nnls")
parser.add_argument("job_id")
args = parser.parse_args()
nr_healthy = int(args.nr_of_healthy_samples)
nr_unhealthy = int(args.nr_unhealthy)
bayes_or_ciber_and_nnls = args.bayes_or_ciber_and_nnls
job_id = int(args.job_id) 
if(bayes_or_ciber_and_nnls == "bayes"):
    number_of_gen_cells_per_cell_state = int(input("please enter a number, that will prepresent the number of genereated cells per cell state, at least 16 or 46 are recommended by the Bayesprism paper"))

print(f"job {job_id} running: inlogged with {nr_healthy} healthy samples, and {nr_unhealthy} unhealthy samples being generated for {bayes_or_ciber_and_nnls}")




def create_rand_matrix():
    current_rand_sig_nr = job_id

    df = pd.read_excel('41590_2017_BFni3693_MOESM10_ESM.xlsx', engine='openpyxl') 
    df.set_index('Majority protein IDs', inplace=True) 
    df_LFQ = df.loc[:, df.columns.str.contains('LFQ.intensity')]
    imputed_df = df.loc[:, ~df.columns.str.contains("Thrombo|Erythro|activ") & df.columns.str.contains("imputed")]
    imputed_df_nonNA = imputed_df.dropna()
    imputed_df_nonNA_sigmatrix04 = imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("04")] 
    imputed_df_nonNA_sigmatrix03 = imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("03")] 
    imputed_df_nonNA_sigmatrix02 = imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("02")] 
    imputed_df_nonNA_sigmatrix01 = imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("01")]

    imputed_df_nonNA_sigmatrix04_inlogged = np.log1p(imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("04")])
    imputed_df_nonNA_sigmatrix03_inlogged = np.log1p(imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("03")])
    imputed_df_nonNA_sigmatrix02_inlogged = np.log1p(imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("02")])
    imputed_df_nonNA_sigmatrix01_inlogged = np.log1p(imputed_df_nonNA.loc[:, imputed_df_nonNA.columns.str.contains("01")])

    n_rows, n_cols = imputed_df_nonNA_sigmatrix01.shape

    # Dirichlet parameters (uniform)
    alpha = [1, 1, 1, 1]
    alpha_real = [1, 1, 1] #this is for creating a sig matrix out of 3 and then testing the left out one, further changes to the code needed if this option is choosen. 

    # Sample Dirichlet weights for each column (shape: [n_cols, 4])
    weights = np.random.dirichlet(alpha, size=n_cols)  # shape: (n_cols, 4)
    weights_real = np.random.dirichlet(alpha_real, size=n_cols)

    # Combine the DataFrames using column-wise Dirichlet weights
    combined_values = (
        weights[:, 0] * imputed_df_nonNA_sigmatrix01.values +
        weights[:, 1] * imputed_df_nonNA_sigmatrix02.values +
        weights[:, 2] * imputed_df_nonNA_sigmatrix03.values +
        weights[:, 3] * imputed_df_nonNA_sigmatrix04.values
    )
    combined_values_inlogged = (
        weights[:, 0] * imputed_df_nonNA_sigmatrix01_inlogged.values +
        weights[:, 1] * imputed_df_nonNA_sigmatrix02_inlogged.values +
        weights[:, 2] * imputed_df_nonNA_sigmatrix03_inlogged.values +
        weights[:, 3] * imputed_df_nonNA_sigmatrix04_inlogged.values
    )
    combined_values_real = (
        weights_real[:, 0] * imputed_df_nonNA_sigmatrix01.values +
        weights_real[:, 1] * imputed_df_nonNA_sigmatrix02.values +
        weights_real[:, 2] * imputed_df_nonNA_sigmatrix03.values
    )
    combined_values_outlogged = np.log1p(combined_values)

    #print(combined_values)
    # Wrap into a DataFrame, preserving original structure
    combined_df = pd.DataFrame(
        combined_values,
        columns=imputed_df_nonNA_sigmatrix04.columns,
        index=imputed_df_nonNA_sigmatrix04.index
    )
    combined_values_inlogged_df = pd.DataFrame(
        combined_values_inlogged, 
        columns=imputed_df_nonNA_sigmatrix04.columns,
        index=imputed_df_nonNA_sigmatrix04.index
    )
    combined_df_outlogged_df = pd.DataFrame(
        combined_values_outlogged,
        columns=imputed_df_nonNA_sigmatrix04.columns,
        index=imputed_df_nonNA_sigmatrix04.index)
    

    # combined_df_real = pd.DataFrame(
    #     combined_values_real, 
    #     columns=imputed_df_nonNA_sigmatrix04.columns,
    #     index=imputed_df_nonNA_sigmatrix04.index
    # )
    combined_df.to_csv(f"imputed_sig_matrix_nonlogged.txt", sep = "\t")
    #combined_df_real.to_csv(f"imputed_sig_matrix_{args}.txt", sep = "\t") #for now....
    combined_df.to_csv(f"sigs/imputed_sig_matrix_{current_rand_sig_nr}.txt", sep = "\t")
    combined_df_outlogged_df.to_csv(f"imputed_sig_matrix_outlogged.txt", sep = "\t")
    combined_values_inlogged_df.to_csv(f"imputed_sig_matrix_inlogged.txt", sep = "\t")



    return
columns_steady_state = [
    "LFQ.intensity_B.memory_01_steady-state",
    "LFQ.intensity_B.memory_02_steady-state",
    "LFQ.intensity_B.memory_03_steady-state",
    "LFQ.intensity_B.memory_04_steady-state",
    "LFQ.intensity_B.naive_01_steady-state",
    "LFQ.intensity_B.naive_02_steady-state",
    "LFQ.intensity_B.naive_03_steady-state",
    "LFQ.intensity_B.naive_04_steady-state",
    "LFQ.intensity_B.plasma_01_steady-state",
    "LFQ.intensity_B.plasma_02_steady-state",
    "LFQ.intensity_B.plasma_03_steady-state",
    "LFQ.intensity_B.plasma_04_steady-state",
    "LFQ.intensity_mTregs_01_steady-state",
    "LFQ.intensity_mTregs_02_steady-state",
    "LFQ.intensity_mTregs_03_steady-state",
    "LFQ.intensity_mTregs_04_steady-state",
    "LFQ.intensity_T4.naive_01_steady-state",
    "LFQ.intensity_T4.naive_02_steady-state",
    "LFQ.intensity_T4.naive_03_steady-state",
    "LFQ.intensity_T4.naive_04_steady-state",
    "LFQ.intensity_nTregs_01_steady-state",
    "LFQ.intensity_nTregs_02_steady-state",
    "LFQ.intensity_nTregs_03_steady-state",
    "LFQ.intensity_nTregs_04_steady-state",
    "LFQ.intensity_T4.CM_01_steady-state",
    "LFQ.intensity_T4.CM_02_steady-state",
    "LFQ.intensity_T4.CM_03_steady-state",
    "LFQ.intensity_T4.CM_04_steady-state",
    "LFQ.intensity_T4.EM_01_steady-state",
    "LFQ.intensity_T4.EM_02_steady-state",
    "LFQ.intensity_T4.EM_03_steady-state",
    "LFQ.intensity_T4.EM_04_steady-state",
    "LFQ.intensity_T4.EMRA_01_steady-state",
    "LFQ.intensity_T4.EMRA_02_steady-state",
    "LFQ.intensity_T4.EMRA_03_steady-state",
    "LFQ.intensity_T4.EMRA_04_steady-state",
    "LFQ.intensity_Th1_01_steady-state",
    "LFQ.intensity_Th1_02_steady-state",
    "LFQ.intensity_Th1_03_steady-state",
    "LFQ.intensity_Th1_04_steady-state",
    "LFQ.intensity_Th17_01_steady-state",
    "LFQ.intensity_Th17_02_steady-state",
    "LFQ.intensity_Th17_03_steady-state",
    "LFQ.intensity_Th17_04_steady-state",
    "LFQ.intensity_Th2_01_steady-state",
    "LFQ.intensity_Th2_02_steady-state",
    "LFQ.intensity_Th2_03_steady-state",
    "LFQ.intensity_Th2_04_steady-state",
    "LFQ.intensity_T8.naive_01_steady-state",
    "LFQ.intensity_T8.naive_02_steady-state",
    "LFQ.intensity_T8.naive_03_steady-state",
    "LFQ.intensity_T8.naive_04_steady-state",
    "LFQ.intensity_T8.CM_01_steady-state",
    "LFQ.intensity_T8.CM_02_steady-state",
    "LFQ.intensity_T8.CM_03_steady-state",
    "LFQ.intensity_T8.CM_04_steady-state",
    "LFQ.intensity_T8.EM_01_steady-state",
    "LFQ.intensity_T8.EM_02_steady-state",
    "LFQ.intensity_T8.EM_03_steady-state",
    "LFQ.intensity_T8.EM_04_steady-state",
    "LFQ.intensity_T8.EMRA_01_steady-state",
    "LFQ.intensity_T8.EMRA_02_steady-state",
    "LFQ.intensity_T8.EMRA_03_steady-state",
    "LFQ.intensity_T8.EMRA_04_steady-state",
    "LFQ.intensity_mDC_01_steady-state",
    "LFQ.intensity_mDC_02_steady-state",
    "LFQ.intensity_mDC_03_steady-state",
    "LFQ.intensity_mDC_04_steady-state",
    "LFQ.intensity_pDC_01_steady-state",
    "LFQ.intensity_pDC_02_steady-state",
    "LFQ.intensity_pDC_03_steady-state",
    "LFQ.intensity_pDC_04_steady-state",
    "LFQ.intensity_Basophil_01_steady-state",
    "LFQ.intensity_Basophil_02_steady-state",
    "LFQ.intensity_Basophil_03_steady-state",
    "LFQ.intensity_Basophil_04_steady-state",
    "LFQ.intensity_Eosinophil_01_steady-state",
    "LFQ.intensity_Eosinophil_02_steady-state",
    "LFQ.intensity_Eosinophil_03_steady-state",
    "LFQ.intensity_Eosinophil_04_steady-state",
    "LFQ.intensity_Neutrophil_01_steady-state",
    "LFQ.intensity_Neutrophil_02_steady-state",
    "LFQ.intensity_Neutrophil_03_steady-state",
    "LFQ.intensity_Neutrophil_04_steady-state",
    "LFQ.intensity_MO.classical_01_steady-state",
    "LFQ.intensity_MO.classical_02_steady-state",
    "LFQ.intensity_MO.classical_03_steady-state",
    "LFQ.intensity_MO.classical_04_steady-state",
    "LFQ.intensity_MO.intermediate_01_steady-state",
    "LFQ.intensity_MO.intermediate_02_steady-state",
    "LFQ.intensity_MO.intermediate_03_steady-state",
    "LFQ.intensity_MO.intermediate_04_steady-state",
    "LFQ.intensity_MO.nonclassical_01_steady-state",
    "LFQ.intensity_MO.nonclassical_02_steady-state",
    "LFQ.intensity_MO.nonclassical_03_steady-state",
    "LFQ.intensity_MO.nonclassical_04_steady-state",
    "LFQ.intensity_NK.bright_01_steady-state",
    "LFQ.intensity_NK.bright_02_steady-state",
    "LFQ.intensity_NK.bright_03_steady-state",
    "LFQ.intensity_NK.bright_04_steady-state",
    "LFQ.intensity_NK.dim_01_steady-state",
    "LFQ.intensity_NK.dim_02_steady-state",
    "LFQ.intensity_NK.dim_03_steady-state",
    "LFQ.intensity_NK.dim_04_steady-state",
    "LFQ.intensity_Erythrocyte_01_steady-state",
    "LFQ.intensity_Erythrocyte_02_steady-state",
    "LFQ.intensity_Erythrocyte_03_steady-state",
    "LFQ.intensity_Erythrocyte_04_steady-state",
    "LFQ.intensity_Thrombocyte_01_steady-state",
    "LFQ.intensity_Thrombocyte_02_steady-state",
    "LFQ.intensity_Thrombocyte_03_steady-state",
    "LFQ.intensity_Thrombocyte_04_steady-state"
]
# Reading input data and first processing into dataframe
df = pd.read_excel('41590_2017_BFni3693_MOESM10_ESM.xlsx', engine='openpyxl') 
df.set_index('Majority protein IDs', inplace=True) 
df_LFQ = df.loc[:, df.columns.str.contains('LFQ.intensity')] #this is the one we will be using
mask_steady = df_LFQ.columns.str.contains('steady')
mask_both = ~df_LFQ.columns.str.contains('imputed')
mask_against_trombo_and_erithro_and_imputed = ~df_LFQ.columns.str.contains('Thrombocyte|Erythrocyte|imputed')
mask_steady_no_thrombo =  mask_against_trombo_and_erithro_and_imputed & mask_steady
df_LFQ_steady_state = df_LFQ.loc[:,mask_steady]
mask = ~df_LFQ_steady_state.columns.str.contains('imputed')
df_both_no_trombo_and_erithro = df_LFQ.loc[:, mask_against_trombo_and_erithro_and_imputed]
df_both_no_trombo_and_erithro_and_imputed = df_LFQ.loc[:, mask_against_trombo_and_erithro_and_imputed] # this one is probs used for sc.dat and so on 
df_LFQ_steady_state = df_LFQ_steady_state.loc[:,mask]
df_LFQ_activated_and_steady_state = df_LFQ.loc[:, mask_both]
df_LFQ_steady_state_no_thrombo = df_LFQ.loc[:, mask_steady_no_thrombo]
imputed_df = df.loc[:, ~df.columns.str.contains("Thrombo|Erythro|activ") & df.columns.str.contains("imputed")]
imputed_df_nonNA = imputed_df.dropna() #some rows were just empty, those were the rows that simply had zero values
# in here we now have all samples, activated and not activated, and some empty lines, how do I get rid of empty lines?
amount_of_that_cell_per_sample = { #from the paper
    'B.memory': 2, 'B.naive': 2.4, 'B.plasma': 0.03, 'Basophil': 1, 'Eosinophil': 0.8,
    'MO.classical': 3.3, 'MO.intermediate': 0.5, 'MO.nonclassical': 1.1, 'NK.bright': 1.6, 'NK.dim': 2.3,
    'Neutrophil': 4.0, 'T4.CM': 3.5, 'T4.EM': 1.7, 'T4.EMRA': 0.2, 'T4.naive': 3, 'T8.CM': 2.1, 'T8.EM': 2.3,
    'T8.EMRA': 1.6, 'T8.naive': 2.5, 'Th1': 2.4, 'Th17': 1.4, 'Th2': 1.5, 'mDC': 0.3,
    'mTregs': 1.1, 'nTregs': 0.4, 'pDC': 0.7
}
dist_myeloid_lymphocytes = pd.read_csv('dist_leukocytes_lymphocytes.csv')
dist_lymphocytes = pd.read_csv('dist_lymphocytes.csv')
#print( dist_lymphocytes.columns)
dist_lymphocytes["Myeloid"] = dist_myeloid_lymphocytes["Monocytes [%] Median"] + dist_myeloid_lymphocytes["PMN [%] Median"]
#dist_lymphocytes["Range Myeloid Bottom"] = dist_myeloid_lymphocytes["Monocytes [%] Range Bottom"]+ dist_myeloid_lymphocytes["PMN [%] Range Bottom"]
#dist_lymphocytes["Range Myeloid Top"] = dist_myeloid_lymphocytes["Monocytes [%] Range Top"] + dist_myeloid_lymphocytes["PMN [%] Range Top"]
dist_lymphocytes.index = dist_lymphocytes["Lymphocyte subset"]
dist_lymphocytes["T"] =  dist_lymphocytes["T helper CD3+4+"] + dist_lymphocytes["T CD3+"] + dist_lymphocytes["T suppressor CD3+8+"]
dist_lymphocytes = dist_lymphocytes.drop(columns=['T helper CD3+4+', 'T CD3+', 'T suppressor CD3+8+'])
dist_lymphocytes = dist_lymphocytes.drop(columns=['Lymphocyte subset', 'Phenotype'])
dist_lymphocytes = dist_lymphocytes.loc[:, ~dist_lymphocytes.columns.str.contains("Range|:")]
dist_lymphocytes = dist_lymphocytes.div(dist_lymphocytes.sum(axis=1), axis=0) *100

def get_generation_samples_sts(cell_type): #returns a dataframe of the 3 samples of one type
    mask = df_LFQ_steady_state.columns.str.contains(cell_type) & ~df_LFQ_steady_state.columns.str.contains('_01_')
    return df_LFQ_steady_state.loc[:, mask]
def get_test_samples_sts(cell_type): #returns a dataframe of the left out test sample of that time
    mask = df_LFQ_steady_state.columns.str.contains(cell_type) & df_LFQ_steady_state.columns.str.contains('_01_')
    return df_LFQ_steady_state.loc[:, mask]
# def get_geneation_samples_sts_imputed_nonNA(cell_type):
#     mask = imputed_df_nonNA.columns.str.contains(cell_type)
#     res = imputed_df_nonNA.loc[:, mask]
#     if (loggedness== "inlogged"):
#         res = np.log1p(res)
    return res #usually should log this
def get_generation_samples_sts_inlogged(cell_type):
    return np.log1p(imputed_df_nonNA.loc[:,imputed_df_nonNA.columns.str.contains(cell_type)])
def get_generation_samples_sts_outlogged_and_nonlogged(cell_type):
    mask = imputed_df_nonNA.columns.str.contains(cell_type)
    return imputed_df_nonNA.loc[:, mask]


cell_type_to_base_mapping = {
    'B.memory': 'B cells',
    'B.naive': 'B cells',
    'B.plasma': 'B cells',
    'Basophil': 'myeloid',
    'Eosinophil': 'myeloid',
    'MO.classical': 'myeloid',
    'MO.intermediate': 'myeloid',
    'MO.nonclassical': 'myeloid',
    'NK.bright': 'NK cells',
    'NK.dim': 'NK cells',
    'Neutrophil': 'myeloid',
    'T4.CM': 'T cells',
    'T4.EM': 'T cells',
    'T4.EMRA': 'T cells',
    'T4.naive': 'T cells',
    'T8.CM': 'T cells',
    'T8.EM': 'T cells',
    'T8.EMRA': 'T cells',
    'T8.naive': 'T cells',
    'Th1': 'T cells',
    'Th17': 'T cells',
    'Th2': 'T cells',
    'mDC': 'myeloid',
    'mTregs': 'T cells',
    'nTregs': 'T cells',
    'pDC': 'myeloid'
}
cell_types = amount_of_that_cell_per_sample.keys()
total_sum = sum(amount_of_that_cell_per_sample.values())
celltype_fraction = {cell_type: amount/total_sum *100 for cell_type, amount in amount_of_that_cell_per_sample.items()}
df_celltype_fraction = pd.read_csv('celltype_fraction.csv')
df_celltype_fraction_percents = df_celltype_fraction.div(df_celltype_fraction.sum(axis=1), axis=0)*100

base_cell_types_mapping  = {
    'B.memory': 'B cells',
    'B.naive': 'B cells',
    'B.plasma': 'B cells',
    'Basophil': 'myeloid',
    'Eosinophil': 'myeloid',
    'MO.classical': 'myeloid',
    'MO.intermediate': 'myeloid',
    'MO.nonclassical': 'myeloid',
    'NK.bright': 'NK cells',
    'NK.dim': 'NK cells',
    'Neutrophil': 'myeloid',
    'T4.CM': 'T cells',
    'T4.EM': 'T cells',
    'T4.EMRA': 'T cells',
    'T4.naive': 'T cells',
    'T8.CM': 'T cells',
    'T8.EM': 'T cells',
    'T8.EMRA': 'T cells',
    'T8.naive': 'T cells',
    'Th1': 'T cells',
    'Th17': 'T cells',
    'Th2': 'T cells',
    'mDC': 'myeloid',
    'mTregs': 'T cells',
    'nTregs': 'T cells',
    'pDC': 'myeloid',
    'B-cell': 'B cells',
    'CD4 memory': 'T cells',
    'CD4 naive': 'T cells',
    'CD8 central memory': 'T cells',
    'CD8 effector memory': 'T cells',
    'CD8 naive': 'T cells',
    'NK CD16+ CD56dim early': 'NK cells',
    'NK CD16+ CD56dim late': 'NK cells',
    'NK CD16- CD56bright': 'NK cells',
    'T-reg': 'T cells',
    'TEMRA': 'T cells',
    'conventional DC2': 'myeloid',
    'dead cell': 'other',
    'dg or MAIT': 'T cells',
    'erithroid': 'other',
    'hspc': 'other',
    'monocyte 1 (CD14+)': 'myeloid',
    'monocyte 2 (CD16+)': 'myeloid',
    'plasma cell': 'B cells',
    'plasmacytoide DC': 'myeloid'
}
adata_celltype_fractions = ad.AnnData(df_celltype_fraction_percents.T, obs= pd.DataFrame(index = df_celltype_fraction_percents.T.index) , var = pd.DataFrame(index = df_celltype_fraction_percents.T.columns))
adata_celltype_fractions.obs['base_cell_type'] = adata_celltype_fractions.obs.index.map(lambda x: base_cell_types_mapping[x])
summed_celltype_fractions = adata_celltype_fractions.to_df().groupby(adata_celltype_fractions.obs['base_cell_type']).sum()

#print("these are the colums" ,dist_lymphocytes.columns)

#now comes the same for coarse type data from tobi
coarse_type_dataframe = pd.read_csv("coarsetype_fraction.csv")
def get_normalized_average_healthy_fraction():  # Returns normalized average fractions for cell types
    result = {}
    for cell_type in dist_lymphocytes.columns:
        # Compute the average fraction across all age groups for the current cell type
        avg_fraction = dist_lymphocytes[cell_type].mean()
        result[cell_type] = avg_fraction
    
    # Normalize the fractions so that their sum equals 1
    total_sum = sum(result.values())
    normalized_result = {cell_type: value / total_sum for cell_type, value in result.items()}
    
    return normalized_result
def get_rand_healthy_frac():#gives a random combination of that celltype fraction #celltypes of the givennn type, it is already coarse
    result =  {}
    for cell_type in dist_lymphocytes.columns:
        #for this specific cell type
        all_params = [dist_lymphocytes.loc[age_group][cell_type] for age_group in dist_lymphocytes.index]
        #11 random numbers adding up to 1
        rand_params = np.random.dirichlet(np.ones(len(all_params)),size=1)
        perccount = [rand_params[0][i]*all_params[i] for i in range(len(all_params))]
        result[cell_type] = np.sum(perccount)
    sum1 = np.sum(list(result.values()))
    result = {cell_type : value / sum1 for cell_type, value in result.items()}
    #print("\n\n this is for the Healthy: \n\n", result)
    return result #returns a dictionary form NK, B, Myeloid, T -> fraction (in that order)
def get_rand_unhealthy_frac(): #same as other funciton just from pandas dataframe, called summed_celltype_fractions, where every column is a sample we want to randomly take the convex combination of
    result =  {}
    for cell_type in summed_celltype_fractions.index:
        #only do this if celltype is not "other"
        if(cell_type != "other"):
            all_params = summed_celltype_fractions.loc[cell_type]
            rand_params = np.random.dirichlet(np.ones(len(all_params)),size=1)
            perccount = [rand_params[0][i]*all_params[i] for i in range(len(all_params))]
            result[cell_type] = np.sum(perccount)
        
    sum1 = np.sum(list(result.values()))
    result = {cell_type : value / sum1 for cell_type, value in result.items()} 
        #reorder such that it makes sense
    consistent_order = ["NK cells", "B cells", "myeloid", "T cells"]
    consistent_reordered_result = {key: result[key] for key in consistent_order}
    #print("\n\n this is for the UNhealthy: \n\n", consistent_reordered_result)
    return consistent_reordered_result #returns a dictionary form B, NK, T, Myeloid Cells -> fraction (in that order)
groups = {
    "B cells": ['B.memory', 'B.naive', 'B.plasma'],
    "T cells": ['T4.CM', 'T4.EM', 'T4.EMRA', 'T4.naive', 'T8.CM', 'T8.EM', 'T8.EMRA', 'T8.naive', 'Th1', 'Th17', 'Th2', 'mTregs', 'nTregs'],
    "Myeloid cells": ['Basophil', 'Eosinophil', 'MO.classical', 'MO.intermediate', 'MO.nonclassical', 'Neutrophil', 'mDC', 'pDC'],
    "NK cells": ['NK.bright', 'NK.dim']
}

# Calculate percentages for each group
percentages = {}
for group_name, cell_types in groups.items():
    # Sum the total for the group
    total = sum(amount_of_that_cell_per_sample[cell_type] for cell_type in cell_types)
    
    # Calculate the percentage for each cell type in the group
    percentages[group_name] = {
        cell_type: (amount_of_that_cell_per_sample[cell_type] / total) * 100 for cell_type in cell_types
    }
def add_to_sample(coarse_group_name, total_sample, coarsetype_sample, healthy, rand_frac): #returns the total sample
    d4array = rand_frac
    if(healthy):
        match coarse_group_name:
            case "B cells":
                total_sample += coarsetype_sample * d4array["B CD19+"]
            case "T cells":
                total_sample += coarsetype_sample * d4array["T"]
            case "Myeloid cells":
                total_sample += coarsetype_sample * d4array["Myeloid"]
            case "NK cells":
                total_sample += coarsetype_sample * d4array["NK CD16.56+3-"]
    if(not healthy):
        match coarse_group_name:
            case "B cells":
                total_sample += coarsetype_sample * d4array["B cells"]
            case "T cells":
                total_sample += coarsetype_sample * d4array["T cells"]
            case "Myeloid cells":
                total_sample += coarsetype_sample * d4array["myeloid"]
            case "NK cells":
                total_sample += coarsetype_sample * d4array["NK cells"]
    return total_sample
all_numpy_arrays = {}
all_imputed_arrays = {}
keyset = {
    'B.memory': 2, 'B.naive': 2.4, 'B.plasma': 0.03, 'Basophil': 1, 'Eosinophil': 0.8,
    'MO.classical': 3.3, 'MO.intermediate': 0.5, 'MO.nonclassical': 1.1, 'NK.bright': 1.6, 'NK.dim': 2.3,
    'Neutrophil': 4.0, 'T4.CM': 3.5, 'T4.EM': 1.7, 'T4.EMRA': 0.2, 'T4.naive': 3, 'T8.CM': 2.1, 'T8.EM': 2.3,
    'T8.EMRA': 1.6, 'T8.naive': 2.5, 'Th1': 2.4, 'Th17': 1.4, 'Th2': 1.5, 'mDC': 0.3,
    'mTregs': 1.1, 'nTregs': 0.4, 'pDC': 0.7}.keys()
for cell_type in keyset:
    selected_columns = get_generation_samples_sts(cell_type + "_")
    selected_columns_imputed = get_generation_samples_sts_inlogged (cell_type + "_")
    numpy_array = selected_columns.to_numpy()
    misterious_numpy_array = selected_columns_imputed
    all_numpy_arrays[cell_type] = numpy_array
    all_imputed_arrays[cell_type] = misterious_numpy_array # here just reformat the dataframe into dictionaries for some reason

def generate_a_sample(): # this is from the non randomized social paper solution
    sample = np.zeros((all_numpy_arrays["B.memory"].shape[0],)) 
    for cell_type in cell_types:
        random_numbers = np.random.dirichlet(np.ones(4), size=1)[0]
        for i in range(4):
            sample += all_numpy_arrays[cell_type][:, i] * random_numbers[i]
    return sample * amount_of_that_cell_per_sample[cell_type]/total_sum 
def generate_coarse_sample(healthy): #we will be only assuming steady state cells at the moment
    sample = np.zeros((all_numpy_arrays["B.memory"].shape[0])) #init sample
    for group_name, group_percentages in percentages.items(): #for each coarse group and their percentages
        basetype_proteins = np.zeros((all_numpy_arrays["B.memory"].shape[0]))
        for cell_type, percentage in group_percentages.items(): #for each subgroup and their percentages
            triple_dataframe = get_generation_samples_sts(cell_type+"_")
            rand_nrs= np.random.dirichlet(np.ones(3), size=1)[0]
            prot_fine_celltype = np.zeros((all_numpy_arrays["B.memory"].shape[0]))
            for i in range(3):
                prot_fine_celltype += triple_dataframe.iloc[:,i] * rand_nrs[i]
            #now we summed up the convex combination of the three samples, what we want is the subtypes to sum up to the coarse type
            basetype_proteins += prot_fine_celltype * percentage
        sample = add_to_sample(group_name, sample, basetype_proteins, healthy)
    return sample
def generate_imputed_sample(healthy):
    sample = np.zeros((all_imputed_arrays['B.naive'].shape[0]))#initializing the sample with the right size
    sample_inlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
    sample_outlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
    sample_nonlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
    #for every coarse cell group
    if healthy: 
        rand_frac = get_rand_healthy_frac() #generates a coarse frac
    else:  
        rand_frac = get_rand_unhealthy_frac()
    for coarse_group, group_percentages in percentages.items():
        coarse_sample = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
        coarse_sample_inlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
        coarse_sample_outlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
        coarse_sample_nonlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
        for fine_celltype, fine_in_coarse_percentage in group_percentages.items():
            # mask = imputed_df_nonNA.columns.str.contains(fine_celltype) & imputed_df_nonNA.columns.str.contains('_04_')
            # fine_type_sample = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
            # fine_type_sample += imputed_df_nonNA.loc[:, mask].iloc[:, 0]
            
            
            triple_df = get_generation_samples_sts_inlogged(fine_celltype+"_")#this name in_fact should be quadruple df
            quadruple_df_inlogged = get_generation_samples_sts_inlogged(fine_celltype+"_")
            quadruple_df_outlogged = get_generation_samples_sts_outlogged_and_nonlogged(fine_celltype+"_")
            quadruple_df_nonlogged = get_generation_samples_sts_outlogged_and_nonlogged(fine_celltype+"_")
            rand_nrs = np.random.dirichlet(np.ones(4), size=1)[0] #change this back to 4
            fine_type_sample = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
            fine_type_sample_inlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
            fine_type_sample_outlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
            fine_type_sample_nonlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
            for i in range(4):
                debugstep = triple_df.iloc[:,i]
                fine_type_sample+=rand_nrs[i]* debugstep
                fine_type_sample_inlogged += rand_nrs[i]*quadruple_df_inlogged.iloc[:,i]
                fine_type_sample_nonlogged += rand_nrs[i] *quadruple_df_nonlogged.iloc[:,i]
                fine_type_sample_outlogged += rand_nrs[i]* quadruple_df_outlogged.iloc[:, i]

        
            coarse_sample+=fine_type_sample* fine_in_coarse_percentage
            coarse_sample_inlogged+=fine_type_sample_inlogged * fine_in_coarse_percentage
            coarse_sample_outlogged+= fine_type_sample_outlogged*fine_in_coarse_percentage
            coarse_sample_nonlogged+=fine_type_sample*fine_in_coarse_percentage


        #figure out what the randomized percentage should be for this celltype
        sample = add_to_sample(coarse_group_name=coarse_group, total_sample=sample, coarsetype_sample=coarse_sample, healthy=healthy, rand_frac=rand_frac) #coarsetype_sample=np.log1p(coarse_sample)
        sample_inlogged = add_to_sample(coarse_group_name=coarse_group, total_sample=sample_inlogged, coarsetype_sample=coarse_sample_inlogged, healthy=healthy, rand_frac=rand_frac)
        sample_nonlogged = add_to_sample(coarse_group_name=coarse_group, total_sample=sample_nonlogged, coarsetype_sample=coarse_sample_nonlogged, healthy= healthy, rand_frac = rand_frac)
        sample_outlogged = add_to_sample(coarse_group_name=coarse_group, total_sample=sample_outlogged, coarsetype_sample=coarse_sample_outlogged, healthy = healthy, rand_frac= rand_frac)
    sample *= 100 #this is now basically sample_nonlogged
    sample_inlogged*=100
    sample_outlogged*=100
    sample_outlogged = np.log1p(sample_outlogged)
    sample_nonlogged*=100
    #sample = np.log1p(sample) #remove this after having tested if it works better
    return (sample, sample_outlogged, sample_inlogged, sample_nonlogged, rand_frac)


def from_coarse_to_fine_fracs(healthy, coarsefracs):
    subtype_fractions = {
        "B.memory": 45.15,
        "B.naive": 54.18,
        "B.plasma": 0.68,
        "T4.CM": 14.77,
        "T4.EM": 7.17,
        "T4.EMRA": 0.84,
        "T4.naive": 12.66,
        "T8.CM": 8.86,
        "T8.EM": 9.70,
        "T8.EMRA": 6.75,
        "T8.naive": 10.55,
        "Th1": 10.13,
        "Th17": 5.91,
        "Th2": 6.33,
        "mTregs": 4.64,
        "nTregs": 1.69,
        "Basophil": 8.55,
        "Eosinophil": 6.84,
        "MO.classical": 28.21,
        "MO.intermediate": 4.27,
        "MO.nonclassical": 9.40,
        "Neutrophil": 34.19,
        "mDC": 2.56,
        "pDC": 5.98,
        "NK.bright": 41.03,
        "NK.dim": 58.97
    }
    
    subtype_to_coarse_mapping_healthy = {
        "B.memory": "B CD19+",
        "B.naive": "B CD19+",
        "B.plasma": "B CD19+",
        "T4.CM": "T",
        "T4.EM": "T",
        "T4.EMRA": "T",
        "T4.naive": "T",
        "T8.CM": "T",
        "T8.EM": "T",
        "T8.EMRA": "T",
        "T8.naive": "T",
        "Th1": "T",
        "Th17": "T",
        "Th2": "T",
        "mTregs": "T",
        "nTregs": "T",
        "Basophil": "Myeloid",
        "Eosinophil": "Myeloid",
        "MO.classical": "Myeloid",
        "MO.intermediate": "Myeloid",
        "MO.nonclassical": "Myeloid",
        "Neutrophil": "Myeloid",
        "mDC": "Myeloid",
        "pDC": "Myeloid",
        "NK.bright": "NK CD16.56+3-",
        "NK.dim": "NK CD16.56+3-"
    }
    subtype_to_coarse_mapping_unhealthy = {
        "B.memory": "B cells",
        "B.naive": "B cells",
        "B.plasma": "B cells",
        "T4.CM": "T cells",
        "T4.EM": "T cells",
        "T4.EMRA": "T cells",
        "T4.naive": "T cells",
        "T8.CM": "T cells",
        "T8.EM": "T cells",
        "T8.EMRA": "T cells",
        "T8.naive": "T cells",
        "Th1": "T cells",
        "Th17": "T cells",
        "Th2": "T cells",
        "mTregs": "T cells",
        "nTregs": "T cells",
        "Basophil": "myeloid",
        "Eosinophil": "myeloid",
        "MO.classical": "myeloid",
        "MO.intermediate": "myeloid",
        "MO.nonclassical": "myeloid",
        "Neutrophil": "myeloid",
        "mDC": "myeloid",
        "pDC": "myeloid",
        "NK.bright": "NK cells",
        "NK.dim": "NK cells"
    }
    if healthy:
        subtype_to_coarse_mapping = subtype_to_coarse_mapping_healthy
    else:   
        subtype_to_coarse_mapping = subtype_to_coarse_mapping_unhealthy

    total_fractions_for_this_sample = {}
    for subtype, fraction in subtype_fractions.items():
        coarse_type = subtype_to_coarse_mapping[subtype]
        total_fractions_for_this_sample[subtype] = (fraction / 100) * coarsefracs[coarse_type]
    return total_fractions_for_this_sample
#mask_sig_matrix = vsned_df_LFQ_steady_state_no_thrombo.columns.str.contains("02_")
mask_sig_matrix_imputed = imputed_df_nonNA.columns.str.contains("04") #lets try 03, now lets try 04 and then 01
imputed_df_nonNA_sigmatrix = imputed_df_nonNA.loc[:, mask_sig_matrix_imputed] #should be be log transformed, maybe



#imputed_df_nonNA_sigmatrix = np.log1p(imputed_df_nonNA_sigmatrix)


imputed_df_nonNA_sigmatrix.to_csv('imputed_sig_matrix.tsv', sep = '\t')
#sig matrix for coarse celltypes
# Define the mapping of coarse cell types to their subcell types and weights
percentages = {
    "B cells": {"B.memory": 45.15, "B.naive": 54.18, "B.plasma": 0.68},
    "T cells": {
        "T4.CM": 14.77, "T4.EM": 7.17, "T4.EMRA": 0.84, "T4.naive": 12.66,
        "T8.CM": 8.86, "T8.EM": 9.70, "T8.EMRA": 6.75, "T8.naive": 10.55,
        "Th1": 10.13, "Th17": 5.91, "Th2": 6.33, "mTregs": 4.64, "nTregs": 1.69
    },
    "Myeloid cells": {
        "Basophil": 8.55, "Eosinophil": 6.84, "MO.classical": 28.21,
        "MO.intermediate": 4.27, "MO.nonclassical": 9.40, "Neutrophil": 34.19,
        "mDC": 2.56, "pDC": 5.98
    },
    "NK cells": {"NK.bright": 41.03, "NK.dim": 58.97}
}

# Normalize percentages to sum to 1 (convert percentages to fractions)
for coarse_type, subtypes in percentages.items():
    total_percentage = sum(subtypes.values())
    percentages[coarse_type] = {subtype: weight / total_percentage for subtype, weight in subtypes.items()}
# Compute the weighted average for the coarse cell types
imputed_df_nonNA_sigmatrix_coarse = pd.DataFrame()
for coarse_type, subtypes in percentages.items():
    # Initialize a column for the coarse cell type
    coarse_type_column = np.zeros(imputed_df_nonNA_sigmatrix.shape[0])
    
    # Iterate over subcell types and their normalized weights
    for subtype, weight in subtypes.items():
        matching_columns = imputed_df_nonNA_sigmatrix.columns[imputed_df_nonNA_sigmatrix.columns.str.contains(subtype)]
        if not matching_columns.empty:
            # Add the weighted contribution of the subcell type
            coarse_type_column += imputed_df_nonNA_sigmatrix[matching_columns[0]] * weight
    
    # Add the computed column to the new DataFrame
    imputed_df_nonNA_sigmatrix_coarse[coarse_type] = coarse_type_column
# Save the coarse signature matrix to a file
imputed_df_nonNA_sigmatrix_coarse.to_csv('imputed_sig_matrix_coarse.tsv', sep='\t')
input_file = "imputed_sig_matrix_coarse.tsv"
output_file = "imputed_sig_matrix_coarse.txt"
with open(input_file, "r") as tsv_file:
    content = tsv_file.read()
with open(output_file, "w") as txt_file:
    txt_file.write(content)

#df_stst_no_th_sig_matrix.to_csv('healthy_sig_matrix.tsv', sep = '\t') #now this is changed to vsned
global_samples_dataframe = pd.DataFrame()
def generate_samples(healthy_sample_nr, unhealthy_sample_nr):
    global global_samples_dataframe
    total_fractions_from_samples = [] # first healthy_sample_nr are the 
    total_coarse_fractions_from_samples = []
    all_samples = []
    all_inlogged_samples = []
    all_outlogged_samples = []
    all_nonlogged_samples = []
    samplenames = []
    for i in range(healthy_sample_nr):
        imputed_sample_1, imputed_sample_outlogged, imputed_sample_inlogged, imputed_sample_nonlogged, coarse_type_fractions1 = generate_imputed_sample(True)
        healthy1= "B CD19+" in coarse_type_fractions1
        total_coarse_fractions_from_samples.append(list(coarse_type_fractions1.values())) #we append the values that come from NK, B, Myeloid, T
        total_fractions_from_samples.append(from_coarse_to_fine_fracs(healthy=healthy1, coarsefracs=coarse_type_fractions1)) 
        samplenames.append("healthy_sample_"+ str(i))
        all_samples.append(imputed_sample_1)
        all_inlogged_samples.append(imputed_sample_inlogged)
        all_outlogged_samples.append(imputed_sample_outlogged)
        all_nonlogged_samples.append(imputed_sample_nonlogged)
    for i in range(unhealthy_sample_nr):
        imputed_sample_1,imputed_sample_outlogged, imputed_sample_inlogged, imputed_sample_nonlogged, coarse_type_fractions1 = generate_imputed_sample(False)
        healthy1= "B CD19+" in coarse_type_fractions1
        total_coarse_fractions_from_samples.append(list(coarse_type_fractions1.values()))
        total_fractions_from_samples.append(from_coarse_to_fine_fracs(healthy=healthy1, coarsefracs=coarse_type_fractions1))
        samplenames.append("unhealthy_sample_" + str(i))
        all_samples.append(imputed_sample_1)
        all_inlogged_samples.append(imputed_sample_inlogged)
        all_outlogged_samples.append(imputed_sample_outlogged)
        all_nonlogged_samples.append(imputed_sample_nonlogged)
    dict_real_fracs = dict(zip(samplenames, total_fractions_from_samples))
    dict_real_coarse_fracs = dict(zip(samplenames, total_coarse_fractions_from_samples))
    #print("this is the real dict coarse fracs ", dict_real_coarse_fracs)
    dict_samples = dict(zip(samplenames, all_samples))
    dict_samples_inlogged = dict(zip(samplenames, all_inlogged_samples))
    dict_samples_outlogged = dict(zip(samplenames, all_outlogged_samples))
    dict_samples_nonlogged = dict(zip(samplenames, all_nonlogged_samples))
    samples_inlogged_dataframe = pd.DataFrame(dict_samples_inlogged)
    samples_outlogged_dataframe = pd.DataFrame(dict_samples_outlogged)
    samples_nonlogged_dataframe = pd.DataFrame(dict_samples_nonlogged)
    samples_dataframe = pd.DataFrame(dict_samples)
    global_samples_dataframe = samples_dataframe
    real_fracs_dataframe = pd.DataFrame(dict_real_fracs)
    real_fracs_coarse_dataframe = pd.DataFrame(dict_real_coarse_fracs)
    samples_dataframe.to_csv('healthy_sample_imputed.tsv', sep = '\t') #not actually healthy but I don't want to have to change names everywhere
    samples_inlogged_dataframe.to_csv('sample_inlogged_imputed.txt', sep = "\t")
    samples_outlogged_dataframe.to_csv('sample_outlogged_imputed.txt', sep = "\t")
    samples_nonlogged_dataframe.to_csv('sample_nonlogged_imputed.txt', sep = "\t")
    real_fracs_dataframe.to_csv('real_fracs.tsv', sep = '\t')
    real_fracs_coarse_dataframe.to_csv('real_coarse_fracs.tsv', sep = '\t')
    # with open(input_file, "r") as tsv_file:
    #     content = tsv_file.read()
    # with open(output_file, "w") as txt_file:
    #     txt_file.write(content)
    print("samples generated :D")
    return 0

def add_generated_columns(df, num_new_columns):
    new_columns = {}

    

#this is for bayesprism inputs
def get_bayes_input(healthy_sample_nr, unhealthy_sample_nr):
    pandas2ri.activate()
    global df_both_no_trombo_and_erithro_and_imputed

    #create a csv for everything
    generate_samples(healthy_sample_nr = healthy_sample_nr, unhealthy_sample_nr = unhealthy_sample_nr)
    #csv for sc.dat 
    #this is the signature matrix
        
    create_rand_matrix()
    randomized_sig_matrix = pd.read_csv(f"imputed_sig_matrix_nonlogged.txt", sep = "\t") #it will just be nonlogged, we cannot easily do all of them.. or can we?
    randomized_sig_matrix = randomized_sig_matrix.drop(columns = ["Majority protein IDs"])
    
    big_sig = imputed_df_nonNA

    with open("unique_mapped_non_NA.txt", "r") as gene_columns:
        lines = [line.strip() for line in gene_columns]
    imputed_df_nonNA_sigmatrix.index = lines
    big_sig.index = lines
    randomized_sig_matrix.index=lines
    imputed_df_nonNA_sigmatrix.to_csv("bayes_sig_matrix.tsv", sep= "\t") #we want imputed...
    randomized_sig_matrix.to_csv("randomized_bayes_sig_matrix.tsv", sep = "\t")

    #we actually want imputed_nonNA for the  sig matrix and also our sig matrix as I understand it will again be simpyl one sample 

    bayes_generated_sig_matrix = imputed_df_nonNA.copy()
    new_columns = {}
    for coarse_group, group_percentages in percentages.items():
        for fine_celltype in group_percentages:
            existing_columns_for_this_celltype = [col for col in bayes_generated_sig_matrix.columns if fine_celltype in col]
            max_patient_nr = max(int(col.split("_")[2]) for col in existing_columns_for_this_celltype)
            for i in range(1, number_of_gen_cells_per_cell_state+1):
                new_patient_nr = max_patient_nr +i
                new_column_name = f"LFQ.intensity.imputed_{fine_celltype}_{new_patient_nr:02d}_steady-state"

                #generating a column
                quadruple_df_nonlogged = get_generation_samples_sts_outlogged_and_nonlogged(fine_celltype+"_")
                fine_type_sample_nonlogged = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
                rand_nrs = np.random.dirichlet(np.ones(4), size=1)[0]
                for i in range(4):
                    debugstep = quadruple_df_nonlogged.iloc[:,i]
                    fine_type_sample_nonlogged+=rand_nrs[i]* debugstep
                new_columns[new_column_name] = fine_type_sample_nonlogged

    for new_col_name, new_col_data in new_columns.items():
        cell_type = new_col_name.split("_")[1]
        existing_columns = [col for col in bayes_generated_sig_matrix.columns if cell_type in col]
        last_existing_column_nr = bayes_generated_sig_matrix.columns.get_loc(existing_columns[-1])
        bayes_generated_sig_matrix.insert(last_existing_column_nr+1, new_col_name, new_col_data)
    print("THIS IS BAYES GEN SIG COLUMNS", bayes_generated_sig_matrix.columns) #this is correct

        #now I should have it into bayes generated
        #next step is to to put it into cell states 


    #csv for cell.state.labels
    # with open("cell_state_labels.tsv", "w") as csv:
    #     for cellsample_name in imputed_df_nonNA.columns.tolist():
    #         csv.write(f"{cellsample_name}\n")
    with open("cell_state_labels.tsv", "w") as csv:
        for cellsample_name in imputed_df_nonNA.columns.tolist():
            csv.write(f"{cellsample_name}\n")
    with open("cell_state_labels_gen.csv", "w") as csv:
        for cellsample_name in bayes_generated_sig_matrix.columns.tolist():
            csv.write(f"{cellsample_name}")
    #csv for bk.dat
    
    bayes_samples_dataframe =  global_samples_dataframe # kinda cool we already have this :))
    bayes_samples_dataframe.index = lines
    bayes_samples_dataframe = bayes_samples_dataframe.T
    bayes_samples_dataframe.to_csv("bayes_samples.tsv", sep='\t') #is tab separated correct even?
    #print(imputed_df_nonNA.columns.tolist()) debugging


    #csv for cell.type.labels
    with open("cell_type_labels.csv", "w") as csv:
        for cellsample_name in imputed_df_nonNA.columns.tolist():
            match = re.search(r"LFQ\.intensity\.imputed_(.*?)_\d+_steady-state", cellsample_name)
            if match:
                csv.write(f'{"_".join(match.groups())}\n')
    with open("cell_type_labels.csv", "r") as f:
        cell_type_labels = f.read().splitlines()
    #csv for cell.type.labels _generated
    with open("cell_type_labels_gen.csv", "w") as csv:
        for cell_ppatient_name in bayes_generated_sig_matrix.columns.tolist():
            match = re.search(r"LFQ\.intensity\.imputed_(.*?)_\d+_steady-state", cell_ppatient_name)
            csv.write(f'{"_".join(match.groups())}\n')
            
    with open("cell_type_labels_gen.csv", "r") as f:
        cell_type_labels_gen = f.read().splitlines()

    r_sc_dat = pandas2ri.py2rpy(big_sig.T) #usually should log this (imputed_df_nonNA_sigmatrix.T)randomized_sig_matrix
    r_sc_dat_gen = pandas2ri.py2rpy(bayes_generated_sig_matrix.T)#same, but more
    r_bk_dat = pandas2ri.py2rpy(bayes_samples_dataframe)
    r_cell_state_labels = ro.StrVector(cell_type_labels)
    r_cell_state_labels_gen = ro.StrVector(cell_type_labels_gen)
    r_cell_type_labels = ro.StrVector(cell_type_labels)
    r_cell_type_labels_gen = ro.StrVector(cell_type_labels_gen)#is this correct? same as cell_state_labels_again?

    # Save the data to an R .rdata file
    ro.r.assign("sc.dat", r_sc_dat)
    ro.r.assign("bk.dat", r_bk_dat)
    ro.r.assign("cell.state.labels", r_cell_state_labels)
    ro.r.assign("cell.type.labels", r_cell_type_labels)
    ro.r('save(bk.dat, cell.state.labels, cell.type.labels, sc.dat, file="myinput2.gbm.rdata")')
    ro.r.assign("sc.dat", r_sc_dat_gen)
    ro.r.assign("bk.dat", r_bk_dat)
    ro.r.assign("cell.state.labels", r_cell_state_labels_gen)
    ro.r.assign("cell.type.labels", r_cell_type_labels_gen)
    ro.r('save(bk.dat, cell.state.labels, cell.type.labels, sc.dat, file="myinput_gen.gbm.rdata")')





    
    return
print("these are healthy fracs", get_normalized_average_healthy_fraction())
if(bayes_or_ciber_and_nnls == "ciber_and_nnls"):
    generate_samples(healthy_sample_nr=nr_healthy, unhealthy_sample_nr=nr_unhealthy)#this is for cibersort inputs 
    create_rand_matrix() # "imputed_sig_matrix_{loggedness}.txt"
    print("rand_sig_matrix_generated")
elif(bayes_or_ciber_and_nnls == "bayes"):
    print("constructing bayes things")
    generate_samples(healthy_sample_nr=nr_healthy, unhealthy_sample_nr=nr_unhealthy)#we need the global dataframe to be set up
    get_bayes_input(nr_healthy, nr_unhealthy) #changed here from 5, 5

else:
    print("not a valid method selection")
    sys.exit()
