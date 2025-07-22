import numpy as np
import pandas as pd
num_samples = 10
sig_sample_nr = "04"
cell_types = [
    "B.memory", "B.naive", "B.plasma",
    "mTregs", "T4.naive", "nTregs", "T4.CM", "T4.EM", "T4.EMRA",
    "Th1", "Th17", "Th2",
    "T8.naive", "T8.CM", "T8.EM", "T8.EMRA",
    "mDC", "pDC", "Basophil", "Eosinophil", "Neutrophil",
    "MO.classical", "MO.intermediate", "MO.nonclassical",
    "NK.bright", "NK.dim"
]
columns = [f"LFQ.intensity.imputed_{ct}_{sig_sample_nr}_steady-state" for ct in cell_types]
data = []
for i in range(num_samples):
    rand_fracs = np.random.dirichlet(np.ones(26), size = 1)[0]
    row = [f"sample_{i+1}"] + list(rand_fracs)
    data.append(row)
df = pd.DataFrame(data, columns= ["Mixture"] + columns)
df.to_csv('RANDOM_GUESSING-Results.txt', sep = "\t", index = False)