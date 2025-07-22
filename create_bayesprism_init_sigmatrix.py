#In this code the signature matrix the way bayesprism uses it, is saved to a file called initial_bayes_ref.txt
#this is useful for the purpouse of estimatinig if the updated signature matrix is somewhat better than this one


import numpy as np
import pandas as pd
alpha = [1, 1, 1, 1]
df = pd.read_excel('41590_2017_BFni3693_MOESM10_ESM.xlsx', engine='openpyxl') 
df.set_index('Majority protein IDs', inplace=True) 
df_LFQ = df.loc[:, df.columns.str.contains('LFQ.intensity')]
imputed_df = df.loc[:, ~df.columns.str.contains("Thrombo|Erythro|activ") & df.columns.str.contains("imputed")]
imputed_df_nonNA = imputed_df.dropna()
# Convert all columns to numeric, coercing non-numeric values to NaN
print(np.random.dirichlet(np.ones(4), size=1)[0])
initial_bayes_ref_matrix = pd.DataFrame()
#Real initial bayesprism matrix: 
for i in range(0, imputed_df_nonNA.shape[1], 4):
    columns_to_sum = imputed_df_nonNA.iloc[:, i:i+4]
    initial_bayes_ref_matrix[columns_to_sum.columns[0]] = columns_to_sum.sum(axis = 1)

print(initial_bayes_ref_matrix.columns)
column_names = [
    'B.memory', 'B.naive', 'B.plasma', 'mTregs','T4.naive', 'nTregs', 'T4.CM', 'T4.EM', 'T4.EMRA', 'Th1', 'Th17',
    'Th2', 'T8.naive',
    'T8.CM', 'T8.EM', 'T8.EMRA', 'mDC', 'pDC',  'Basophil', 'Eosinophil','Neutrophil',  'MO.classical',
    'MO.intermediate', 'MO.nonclassical', 
    'NK.bright', 'NK.dim'
]
initial_bayes_ref_matrix.columns = column_names[:initial_bayes_ref_matrix.shape[1]]
initial_bayes_ref_matrix.to_csv("initial_bayes_ref_matrix.txt", sep = "\t")