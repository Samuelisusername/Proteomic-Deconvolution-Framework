#imports
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
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
df_both_no_trombo_and_erithro_and_imputed = df_LFQ.loc[:, mask_against_trombo_and_erithro_and_imputed]
df_LFQ_steady_state = df_LFQ_steady_state.loc[:,mask]
df_LFQ_activated_and_steady_state = df_LFQ.loc[:, mask_both]
df_LFQ_steady_state_no_thrombo = df_LFQ.loc[:, mask_steady_no_thrombo]
imputed_df = df.loc[:, ~df.columns.str.contains("Thrombo|Erythro|activ") & df.columns.str.contains("imputed")]
imputed_df_nonNA = imputed_df.dropna() #some rows were just empty, those were the rows that simply had zero values
imputed_df_nonNA.to_excel("imputed_only_nonNA.xlsx")# in here we now have all samples, activated and not activated, and some empty lines, how do I get rid of empty lines?
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
def get_geneation_samples_sts_imputed_nonNA(cell_type):
    mask = imputed_df_nonNA.columns.str.contains(cell_type) & ~imputed_df_nonNA.columns.str.contains('_01_')
    return np.log1p(imputed_df_nonNA.loc[:, mask]) # this line was changed


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
    print("\n\n this is for the Healthy: \n\n", result)
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
    print("\n\n this is for the UNhealthy: \n\n", consistent_reordered_result)
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
    selected_columns_imputed = get_geneation_samples_sts_imputed_nonNA(cell_type + "_")
    numpy_array = selected_columns.to_numpy()
    misterious_numpy_array = selected_columns_imputed
    all_numpy_arrays[cell_type] = numpy_array
    all_imputed_arrays[cell_type] = misterious_numpy_array # here just reformat the dataframe into dictionaries for some reason

def generate_a_sample(): # this is from the non randomized social paper solution
    sample = np.zeros((all_numpy_arrays["B.memory"].shape[0],)) 
    for cell_type in cell_types:
        random_numbers = np.random.dirichlet(np.ones(3), size=1)[0]
        for i in range(3):
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
    #for every coarse cell group
    if healthy: 
        rand_frac = get_rand_healthy_frac() #generates a coarse frac
    else:  
        rand_frac = get_rand_unhealthy_frac()
    for coarse_group, group_percentages in percentages.items():
        coarse_sample = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
        for fine_celltype, fine_in_coarse_percentage in group_percentages.items():
            triple_df = get_geneation_samples_sts_imputed_nonNA(fine_celltype+"_")
            rand_nrs = np.random.dirichlet(np.ones(3), size=1)[0]
            fine_type_sample = np.zeros((all_imputed_arrays['B.naive'].shape[0]))
            for i in range(3):
                debugstep = triple_df.iloc[:,i]
                fine_type_sample+=rand_nrs[i]* debugstep
            coarse_sample+=fine_type_sample* fine_in_coarse_percentage
        #figure out what the randomized percentage should be for this celltype
        sample = add_to_sample(coarse_group_name=coarse_group, total_sample=sample, coarsetype_sample=coarse_sample, healthy=healthy, rand_frac=rand_frac)
    sample *= 100 
    #sample = np.log1p(sample) #remove this after having tested if it works better
    return (sample, rand_frac)






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
imputed_df_nonNA_sigmatrix = np.log1p(imputed_df_nonNA.loc[:, mask_sig_matrix_imputed])



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
if(imputed_df_nonNA_sigmatrix.index.equals(imputed_df_nonNA_sigmatrix_coarse.index)):
    print("they are the same!!!")
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
    samplenames = []
    for i in range(healthy_sample_nr):
        imputed_sample_1, coarse_type_fractions1 = generate_imputed_sample(True)
        healthy1= "B CD19+" in coarse_type_fractions1
        total_coarse_fractions_from_samples.append(list(coarse_type_fractions1.values())) #we append the values that come from NK, B, Myeloid, T
        total_fractions_from_samples.append(from_coarse_to_fine_fracs(healthy=healthy1, coarsefracs=coarse_type_fractions1)) 
        samplenames.append("healthy_sample_"+ str(i))
        all_samples.append(imputed_sample_1)
    for i in range(unhealthy_sample_nr):
        imputed_sample_1, coarse_type_fractions1 = generate_imputed_sample(False)
        healthy1= "B CD19+" in coarse_type_fractions1
        total_coarse_fractions_from_samples.append(list(coarse_type_fractions1.values()))
        total_fractions_from_samples.append(from_coarse_to_fine_fracs(healthy=healthy1, coarsefracs=coarse_type_fractions1))
        samplenames.append("unhealthy_sample_" + str(i))
        all_samples.append(imputed_sample_1)
    dict_real_fracs = dict(zip(samplenames, total_fractions_from_samples))
    dict_real_coarse_fracs = dict(zip(samplenames, total_coarse_fractions_from_samples))
    #print("this is the real dict coarse fracs ", dict_real_coarse_fracs)
    dict_samples = dict(zip(samplenames, all_samples))
    samples_dataframe = pd.DataFrame(dict_samples)
    global_samples_dataframe = samples_dataframe
    real_fracs_dataframe = pd.DataFrame(dict_real_fracs)
    real_fracs_coarse_dataframe = pd.DataFrame(dict_real_coarse_fracs)
    if samples_dataframe.index.equals(imputed_df_nonNA_sigmatrix_coarse.index):
        print("those are the same tooo")
    samples_dataframe.to_csv('healthy_sample_imputed.tsv', sep = '\t') #not actually healthy but I don't want to have to change names everywhere
    real_fracs_dataframe.to_csv('real_fracs.tsv', sep = '\t')
    real_fracs_coarse_dataframe.to_csv('real_coarse_fracs.tsv', sep = '\t')
    input_file = "healthy_sample_imputed.tsv"
    output_file = "healthy_sample_imputed.txt"
    with open(input_file, "r") as tsv_file:
        content = tsv_file.read()
    with open(output_file, "w") as txt_file:
        txt_file.write(content)
    input_file = "imputed_sig_matrix.tsv"
    output_file = "imputed_sig_matrix.txt"
    with open(input_file, "r") as tsv_file:
        content = tsv_file.read()
    with open(output_file, "w") as txt_file:
        txt_file.write(content)
    print("it worked :D")
    return 0
def generate_coarse_samples(healthy_sample_nr, unhealthy_sample_nr):
    global global_samples_dataframe #here we will store the samples

#this is for bayesprism inputs
def get_bayes_input(healthy_sample_nr, unhealthy_sample_nr):
    global df_both_no_trombo_and_erithro_and_imputed

    #create a csv for everything
    generate_samples(healthy_sample_nr = healthy_sample_nr, unhealthy_sample_nr = unhealthy_sample_nr)
    #csv for sc.dat 
    #this is the signature matrix
    df_both_no_trombo_and_erithro_and_imputed.to_csv("bayes_sig_matrix.tsv", sep= "\t")
    #csv for cell.state.labels
    with open("cell_state_labels.tsv", "w") as csv:
        for cellsample_name in df_both_no_trombo_and_erithro_and_imputed.columns.tolist():
            csv.write(f"{cellsample_name}\n")
    #csv for bk.dat
    bayes_samples_dataframe =  global_samples_dataframe.T # kinda cool we already have this :))
    bayes_samples_dataframe.to_csv("bayes_samples.tsv", sep='\t') #is tab separated correct even?
        
    #csv for cell.type.labels
    with open("cell_type_labels.csv", "w") as csv:
        for cellsample_name in df_both_no_trombo_and_erithro_and_imputed.columns.tolist():
            match = re.search(r"LFQ\.intensity_(.*?)_.*?_(steady-state|activated)", cellsample_name)
            if match:
                csv.write(f'{"_".join(match.groups())}\n')
    return

generate_samples(healthy_sample_nr=5, unhealthy_sample_nr=5)#this is for cibersort inputs
#get_bayes_input(3, 3)
