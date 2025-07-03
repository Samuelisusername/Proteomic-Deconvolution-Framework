healthy = 5
unhealthy = 5
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
data = []
healthy_coarse_fracs = {"B cells":9.125807636613832,"T cells":64.09516066434331, "Myeloid cells":22.551629396455482, "NK cells": 4.227402302587374}
healthy_fine_fracs = [[coarse_frac* fine_frac /100 for fine_frac in percentages[celltype].values()] for (celltype, coarse_frac) in healthy_coarse_fracs.items()]
biglist = []
for list in healthy_fine_fracs:
    for elem in list:
        biglist.append(elem)
print(biglist)