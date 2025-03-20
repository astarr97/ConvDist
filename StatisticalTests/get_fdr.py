import os
import pandas as pd
import numpy as np

folder = "Aquatic_ConvNearest_Results/Combined_Files/"
midfix = "ConvDiv_Combined"

for prefix in ["NC_GOBP", "NC_HPO", "NC_PerGene"]:
    df = pd.DataFrame()
    for file in os.listdir(folder):
        if prefix in file and midfix in file and "Permute" in file:
            #print(file)
            v = pd.read_csv(folder + file)
            #print(v)
            v = v[["Combined p-value all comps Conv", "Combined p-value all comps Div", "ConvDiv Combined p-value"]]
            df = pd.concat([v, df])
        elif prefix in file and midfix in file:
            v_real = pd.read_csv(folder + file)
            real_file = file
    df = df.sort_values("ConvDiv Combined p-value")
    
    tot = df.shape[0]
    num_tests = v_real.shape[0]
    fdr_perm_ConvDiv = []
    fdr_perm_Conv = []
    fdr_perm_Div = []
    
    c = 1
    for index, row in v_real.iterrows():
        fdr_perm_ConvDiv.append(np.min([df[df["ConvDiv Combined p-value"] <= row["ConvDiv Combined p-value"]].shape[0]*num_tests/tot/c, 1]))
        fdr_perm_Div.append(np.min([df[df["Combined p-value all comps Div"] <= row["Combined p-value all comps Div"]].shape[0]*num_tests/tot/c, 1]))
        fdr_perm_Conv.append(np.min([df[df["Combined p-value all comps Conv"] <= row["Combined p-value all comps Conv"]].shape[0]*num_tests/tot/c, 1]))
    
        c += 1
    v_real["ConvDiv Combined FDR Perm"] = fdr_perm_ConvDiv
    v_real["Combined FDR all comps Div Perm"] = fdr_perm_Div
    v_real["Combined FDR all comps Conv Perm"] = fdr_perm_Conv

    v_real.to_csv(real_file.replace(".csv", "_PermFDR.csv"), sep = ",", index = False)