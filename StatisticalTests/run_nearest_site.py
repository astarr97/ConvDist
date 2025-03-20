#import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd
import numpy as np
#import copy
#import seaborn as sns
from scipy.stats import mannwhitneyu as mwu
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import wilcoxon
from scipy.optimize import curve_fit
from scipy.stats import fisher_exact
import os
from scipy.stats import combine_pvalues
from scipy.stats import spearmanr,pearsonr
from collections import Counter
import sys

#Need to define spec_sup, permute (bool), dist_cut, species_dif_cut, file_list, folder_name

hfont = {'fontname':'Arial'}
#plt.rcParams["font.family"] = "Arial"

#Code borrowed heavily from here: https://stackoverflow.com/questions/62375034/find-non-overlapping-area-between-two-kde-plots
#plt.rcParams.update({"text.usetex": False})


file_list = sys.argv[1]
file_list = file_list.split(",")
folder_name = sys.argv[2]
prefix_run = sys.argv[3]
spec_sup = float(sys.argv[4])
dist_cut = float(sys.argv[5])
species_dif_cut = float(sys.argv[6])
permute = int(sys.argv[7])
print(permute)
if permute:
    seed = int(sys.argv[8])

out_dir = folder_name + "_ConvNearest_Results"


def quantile_normalize(v):
    v.index = v["Chrom"] + ":" + v["Pos2"].astype(str)
    
    #Make dataframe to normalize
    df = v[["PhyloP", "Mean PhyloP"]].copy().sort_values("PhyloP")
    df2 = df[["Mean PhyloP"]].sort_values("Mean PhyloP")
    
    df["Mean PhyloP ToNorm"] = list(df2["Mean PhyloP"])
    df.index = df.index + "-" + df2.index
    
    df = df[["PhyloP", "Mean PhyloP ToNorm"]].copy()

    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    df.columns = ["Corrected Quant PhyloP", "Corrected Quant Mean PhyloP"]
    df1 = df[["Corrected Quant PhyloP"]].copy()
    df2 = df[["Corrected Quant Mean PhyloP"]].copy()
    df1.index = [x.split("-")[0] for x in df1.index]
    df2.index = [x.split("-")[1] for x in df2.index]
    
    v_to_ret = v.join(df1)
    v_to_ret = v_to_ret.join(df2)
    return v_to_ret

#Functions to permute the conv/div and nearest sites to ensure the test is well-calibrated
def pni(v, seed):
    np.random.seed(seed)
    new_phylop = []
    new_mean_phylop = []
    for index, row in v.iterrows():
        if np.random.choice([0, 1]):
            new_phylop.append(row["Corrected Quant PhyloP"])
            new_mean_phylop.append(row["Corrected Quant Mean PhyloP"])
        else:
            new_mean_phylop.append(row["Corrected Quant PhyloP"])
            new_phylop.append(row["Corrected Quant Mean PhyloP"])
    v["Corrected Quant PhyloP"] = new_phylop
    v["Corrected Quant Mean PhyloP"] = new_mean_phylop
    return v

def permute_nearest(v_new_cds, v_conv_cds, v_div_cds, v_new_nc, v_conv_nc, v_div_nc, seed=6):
    return pni(v_new_cds, seed), pni(v_conv_cds, seed), pni(v_div_cds, seed), pni(v_new_nc, seed), pni(v_conv_nc, seed), pni(v_div_nc, seed)

def read_ConvNearest_file(file, spec_sup=100, dist_cut=1000, species_dif_cut=10):

    v = pd.read_csv(file, sep = "\t", header = None)
    
    v = v.drop_duplicates([0, 2])
    v["Position1"] = v[0] + ":" + v[2].astype(str)
    v["Position2"] = v[26] + ":" + v[28].astype(str)
    v = v[~v["Position1"].isin(v["Position2"])]
    
    #May need to change this as it may error in the future when there are no "." in the column!
    v = v[v[4] != "."]
    v[3] = v[3].astype(float)
    v = v[v[3] != "."]
    v = v[v[6] != "."]
    v[6] = v[6].astype(int)
    v[5] = v[5].astype(int)

    #21 is PhyloP 1
    v = v[v[21] != "."]
    v[21] = v[21].astype(float)
    v = v[v[22] != "."]
    v = v[v[24] != "."]
    v[24] = v[24].astype(int)
    v[23] = v[23].astype(int)
    v[25] = v[25].astype(int)

    #32 is PhyloP 2
    v = v[v[32] != "."]
    v[32] = v[32].astype(float)
    v = v[v[33] != "."]
    v = v[v[35] != "."]
    v[35] = v[35].astype(int)
    v[34] = v[34].astype(int)
    v[36] = v[36].astype(int)
    
    #After we converted to the proper format
    #We go through restricting to only sites that are AncestralSame (row[19] == row[20] guarantees this)
    #We additionally require that the distance to the convergent/divergent site be < dist_cut 1000 bases is the default (row[25] < dist_cut)
    #And require that the species support be within species_dif_cut (default 10) of the convergent/divergent site
    out = []
    first_seen = []
    second_seen = []
    for index, row in v.iterrows():
        first_good = False
        second_good = False
        if row[14] == "AncestralSame":
            if "Rel" in file:
                to_compare_out1 = 18
                if "Loxodonta_africana" == file.split(".")[0] or "Trichechus_manatus" == file.split(".")[0]:
                    to_compare_out1 = 19
                to_compare_out2 = 29
                #Specifically deal with the fucking manatee clade which does not have a decent outgroup!!!
                if "Loxodonta_africana" == file.split(".")[1] or "Trichechus_manatus" == file.split(".")[1]:
                    to_compare_out2 = 30
            else:
                to_compare_out1 = 19
                to_compare_out2 = 30
            if row[to_compare_out1] == row[20] and row[25] < dist_cut and abs(row[24] - row[6]) < species_dif_cut:
                #This block makes sure the nearest site is in a CDS if the main site is or is not in a CDS if the nearest site is not in a CDS
                if row[5] == 0:
                    if row[23] == 0:
                        first_good = True
                else:
                    if row[23] != 0:
                        first_good = True
                
            if row[to_compare_out2] == row[31] and row[36] < dist_cut and abs(row[35] - row[6]) < species_dif_cut:
                #This block makes sure the nearest site is in a CDS if the main site is or is not in a CDS if the nearest site is not in a CDS
                if row[5] == 0:
                    if row[34] == 0:
                        second_good = True
                else:
                    if row[34] != 0:
                        second_good = True
            #Average if both sites are good, take one or the other if not
            if first_good and second_good:
                second_seen.append(row[26] + ":" + str(row[28]))
                first_seen.append(row[15] + ":" + str(row[17]))
                if row[6] > spec_sup and row[24] > spec_sup and row[35] > spec_sup:
                    out.append([row[0], row[2], row[3], row[4], row[5], row[6], row[21], row[32], np.mean([row[21], row[32]]), row[13]])
            elif first_good:
                second_seen.append(row[15] + ":" + str(row[17]))
                if row[6] > spec_sup and row[24] > spec_sup:
                    out.append([row[0], row[2], row[3], row[4], row[5], row[6], row[21], "NA", row[21], row[13]])
            elif second_good:
                first_seen.append(row[26] + ":" + str(row[28]))
                if row[6] > spec_sup and row[35] > spec_sup:
                    out.append([row[0], row[2], row[3], row[4], row[5], row[6], "NA", row[32], row[32], row[13]])

    v_new = pd.DataFrame(out)
    print(v_new)
    print(file)
    print("Prints out how many nearest sites are used >= 2 times, if this is a high percentage they may be worth removing")
    print(v_new)
    c = Counter(first_seen)
    doubles = 0
    singles = 0
    for key in c.keys():
        if c[key] > 1:
            doubles += 1
        else:
            singles += 1
    print("First species: ", "Proportion doubles: ", doubles*2/singles, "Num doubles: ", doubles, "Num singles: ", singles)
    
    c = Counter(second_seen)
    doubles = 0
    singles = 0
    for key in c.keys():
        if c[key] > 1:
            doubles += 1
        else:
            singles += 1
    print("Second species: ", "Proportion doubles: ", doubles*2/singles, "Num doubles: ", doubles, "Num singles: ", singles)

    
    v_new.columns = ["Chrom", "Pos2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "PhyloP S1", "PhyloP S2", "Mean PhyloP", "ConvOrDiv"]
    v_new = v_new.drop_duplicates(["Chrom", "Pos2"])
    
    #Split into convergent/divergent 
    v_conv = v_new[v_new["ConvOrDiv"].isin(["Convergent"])].copy()
    v_div = v_new[v_new["ConvOrDiv"].isin(["Divergent"])].copy()
    
    print("Number convergent: " + str(len(v_conv.index)))
    print("Number divergent: " + str(len(v_div.index)))
    
    #Split into cds/nc
    v_conv_cds = v_conv[v_conv["NearestDist"] == 0]
    v_div_cds = v_div[v_div["NearestDist"] == 0]
    v_new_cds = v_new[v_new["NearestDist"] == 0]
    
    v_conv_nc = v_conv[v_conv["NearestDist"] != 0]
    v_div_nc = v_div[v_div["NearestDist"] != 0]
    v_new_nc = v_new[v_new["NearestDist"] != 0]
    
    #Quantile normalize
    v_conv_cds["Corrected"] = v_conv_cds["PhyloP"] - (np.mean(v_conv_cds["PhyloP"]) - np.mean(v_conv_cds["Mean PhyloP"]))
    v_div_cds["Corrected"] = v_div_cds["PhyloP"] - (np.mean(v_div_cds["PhyloP"]) - np.mean(v_div_cds["Mean PhyloP"]))
    v_conv_nc["Corrected"] = v_conv_nc["PhyloP"] - (np.mean(v_conv_nc["PhyloP"]) - np.mean(v_conv_nc["Mean PhyloP"]))
    v_div_nc["Corrected"] = v_div_nc["PhyloP"] - (np.mean(v_div_nc["PhyloP"]) - np.mean(v_div_nc["Mean PhyloP"]))

    v_conv_cds = quantile_normalize(v_conv_cds)
    v_conv_nc = quantile_normalize(v_conv_nc)
    v_div_cds = quantile_normalize(v_div_cds)
    v_div_nc = quantile_normalize(v_div_nc)

    v_new_cds["Corrected"] = v_new_cds["PhyloP"] - (np.mean(v_new_cds["PhyloP"]) - np.mean(v_new_cds["Mean PhyloP"]))
    v_new_cds = quantile_normalize(v_new_cds)
    v_new_nc["Corrected"] = v_new_nc["PhyloP"] - (np.mean(v_new_nc["PhyloP"]) - np.mean(v_new_nc["Mean PhyloP"]))
    v_new_nc = quantile_normalize(v_new_nc)
    
    
    return v_new_cds, v_conv_cds, v_div_cds, v_new_nc, v_conv_nc, v_div_nc

#Define function to search for convergent evolution at gene set level
#Default min_var_cut is 25 for gene sets
def run_gene_set_test(v, gene_set, min_var_cut=25, stat_test="Corrected Quant", stat_print="PhyloP"):
    out = []
    out_g0 = []
    if type(gene_set) is dict:
        for key in gene_set.keys():
            v2 = v[v["NearestGene"].isin(gene_set[key])].copy()
            if len(v2.index) > min_var_cut:
                try:
                    out.append([key, v2.shape[0], np.median(v2["PhyloP"]), np.median(v2["Mean " + stat_print]), np.mean(np.array(v2[stat_test + " " + stat_print]) - np.array(v2[stat_test + " Mean " + stat_print])), wilcoxon(v2[stat_test + " " + stat_print], v2[stat_test + " Mean " + stat_print], alternative = "greater")[1]])
                except:
                    print(key)
    elif type(gene_set) is list:
        gene_set = np.unique(v["NearestGene"])
        for gene in gene_set:
            v2 = v[v["NearestGene"].isin([gene])].copy()
            if len(v2.index) > min_var_cut:
                try:
                    out.append([gene, v2.shape[0], np.median(v2["PhyloP"]), np.median(v2["Mean " + stat_print]), np.mean(np.array(v2[stat_test + " " + stat_print]) - np.array(v2[stat_test + " Mean " + stat_print])), wilcoxon(v2[stat_test + " " + stat_print], v2[stat_test + " Mean " + stat_print], alternative = "greater")[1]])
                except:
                    print(gene)


    df = pd.DataFrame(out)
    if df.shape[0]:
        df["FDR"] = fdrcorrection(df[5])[1]
        df = df.sort_values(5)
        df.columns = ["Term", "Number Sites", "Median " + stat_print + " Convergent", "Median " + stat_print + " Nearest", "Paired Mean Difference " + stat_test, "Wilcoxon p-value " + stat_test, "Wilcoxon FDR"]
        
    return df


#Read in gene sets
gobp = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Ontologies/GOBP_AccelEvol_Input.txt", sep= "\t")
d_BP = {}

for index, row in gobp.iterrows():
    d_BP[row["Term"]] = row["Genes"].split(";")
    
hpo = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Ontologies/HPO_AccelEvol_Input.txt", sep= "\t")
d_HPO = {}

for index, row in hpo.iterrows():
    d_HPO[row["Term"]] = row["Genes"].split(";")
    
gene_org = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Ontologies/GeneOrganizer_AccelEvol_Input.txt", sep= "\t")
d_GeneOrg = {}

for index, row in gene_org.iterrows():
    d_GeneOrg[row["Term"]] = row["Genes"].split(";")
    

#To combine p-values, we need to make sure the same sites do not appear in different comparisons
#To deal with this, we remove all sites that have already been tested
if out_dir not in os.listdir():
    os.mkdir(out_dir)
remove_sites = []
for i in range(len(file_list)):
    file = file_list[i]
    v_new_cds, v_conv_cds, v_div_cds, v_new_nc, v_conv_nc, v_div_nc = read_ConvNearest_file(file, species_dif_cut = 10)
    if permute:
        print("Permuting div or conv with nearest")
        v_new_cds, v_conv_cds, v_div_cds, v_new_nc, v_conv_nc, v_div_nc = permute_nearest(v_new_cds, v_conv_cds, v_div_cds, v_new_nc, v_conv_nc, v_div_nc, seed=seed)
    #To run per gene we can just pass the empty list
    genes_cds = []
    genes_nc = []
    
    #Run the tests on HPO, GOBP, and per gene
    map_prefix = {"NC_GOBP_":d_BP, "NC_HPO_":d_HPO, "NC_PerGene_":genes_nc, "CDS_GOBP_":d_BP, "CDS_HPO_":d_HPO, "CDS_PerGene_":genes_nc}
    for key in map_prefix.keys():
        if "CDS" == key[0:3]:
            if map_prefix[key] == []:
                min_var_cut = 5
            else:
                min_var_cut = 10
        elif "NC" == key[0:2]:
            if map_prefix[key] == []:
                min_var_cut = 10
            else:
                min_var_cut = 25
        if "NC" == key[0:2]:
            
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".ConvNearest.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_conv_nc, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".DivNearest.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_div_nc, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            #Comment out for now as the current plan is to combine p-values for Conv and Div sites if desired
            #out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".AllNearest.Prelim.csv"
            #if permute:
            #    out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            #df = run_gene_set_test(v_new_nc, map_prefix[key])
            #if df.shape[0]:
                #df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
        elif "CDS" == key[0:3]:
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".ConvNearest.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_conv_cds, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".DivNearest.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_div_cds, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            #Comment out for now as the current plan is to combine p-values for Conv and Div sites if desired
            #out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".AllNearest.Prelim.csv"
            #if permute:
            #    out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            #df = run_gene_set_test(v_new_cds, map_prefix[key])
            #if df.shape[0]:
                #df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
    
    v_conv_cds["Chrom:Pos2"] = v_conv_cds["Chrom"] + ":" + v_conv_cds["Pos2"].astype(str)
    v_conv_cds = v_conv_cds[~v_conv_cds["Chrom:Pos2"].isin(remove_sites)]
    remove_sites = np.unique(list(remove_sites) + list(v_conv_cds["Chrom:Pos2"]))
    v_conv_cds = v_conv_cds.drop("Chrom:Pos2", axis = 1)
    
    print(len(v_conv_cds.index))
    v_div_cds["Chrom:Pos2"] = v_div_cds["Chrom"] + ":" + v_div_cds["Pos2"].astype(str)
    v_div_cds = v_div_cds[~v_div_cds["Chrom:Pos2"].isin(remove_sites)]
    remove_sites = np.unique(list(remove_sites) + list(v_div_cds["Chrom:Pos2"]))
    v_div_cds = v_div_cds.drop("Chrom:Pos2", axis = 1)
    
    v_conv_nc["Chrom:Pos2"] = v_conv_nc["Chrom"] + ":" + v_conv_nc["Pos2"].astype(str)
    v_conv_nc = v_conv_nc[~v_conv_nc["Chrom:Pos2"].isin(remove_sites)]
    remove_sites = np.unique(list(remove_sites) + list(v_conv_nc["Chrom:Pos2"]))
    v_conv_nc = v_conv_nc.drop("Chrom:Pos2", axis = 1)

    v_div_nc["Chrom:Pos2"] = v_div_nc["Chrom"] + ":" + v_div_nc["Pos2"].astype(str)
    v_div_nc = v_div_nc[~v_div_nc["Chrom:Pos2"].isin(remove_sites)]
    remove_sites = np.unique(list(remove_sites) + list(v_div_nc["Chrom:Pos2"]))
    v_div_nc = v_div_nc.drop("Chrom:Pos2", axis = 1)
    
    #Again, commented out for now
    #v_new_cds["Chrom:Pos2"] = v_new_cds["Chrom"] + ":" + v_new_cds["Pos2"].astype(str)
    #v_new_cds = v_new_cds[~v_new_cds["Chrom:Pos2"].isin(remove_sites)]
    #remove_sites = np.unique(list(remove_sites) + list(v_new_cds["Chrom:Pos2"]))
    #v_new_cds = v_new_cds.drop("Chrom:Pos2", axis = 1)

    #v_new_nc["Chrom:Pos2"] = v_new_nc["Chrom"] + ":" + v_new_nc["Pos2"].astype(str)
    #v_new_nc = v_new_nc[~v_new_nc["Chrom:Pos2"].isin(remove_sites)]
    #remove_sites = np.unique(list(remove_sites) + list(v_new_nc["Chrom:Pos2"]))
    #v_new_nc = v_new_nc.drop("Chrom:Pos2", axis = 1)
    
    for key in map_prefix.keys():
        if "CDS" == key[0:3]:
            if map_prefix[key] == []:
                min_var_cut = 5
            else:
                min_var_cut = 10
        elif "NC" == key[0:2]:
            if map_prefix[key] == []:
                min_var_cut = 10
            else:
                min_var_cut = 25
        if "NC" == key[0:2]:
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".ConvNearest.ForCombine.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_conv_nc, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".DivNearest.ForCombine.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_div_nc, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            #Comment out for now as the current plan is to combine p-values for Conv and Div sites if desired
            #out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".AllNearest.ForCombine.Prelim.csv"
            #if permute:
            #    out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            #df = run_gene_set_test(v_new_nc, map_prefix[key], min_var_cut = min_var_cut)
            #if df.shape[0]:
                #df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
        elif "CDS" == key[0:3]:
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".ConvNearest.ForCombine.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_conv_cds, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".DivNearest.ForCombine.Prelim.csv"
            if permute:
                out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            df = run_gene_set_test(v_div_cds, map_prefix[key], min_var_cut = min_var_cut)
            if df.shape[0]:
                df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
            
            #Comment out for now as the current plan is to combine p-values for Conv and Div sites if desired
            #out_suffix = file.split(".")[0] + "." + file.split(".")[1] + ".AllNearest.ForCombine.Prelim.csv"
            #if permute:
            #    out_suffix = out_suffix.replace(".csv", ".Permute.Seed" + str(seed) + ".csv")
            #df = run_gene_set_test(v_new_cds, map_prefix[key], min_var_cut = min_var_cut)
            #if df.shape[0]:
                #df.sort_values("Wilcoxon p-value Corrected Quant").to_csv(out_dir + "/" + key + out_suffix, index = False)
    
    
    print(len(remove_sites))
   
if "Combined_Files" not in os.listdir(out_dir):
    os.mkdir(out_dir + "/Combined_Files")

#Combine p-values for the convergent sites
out = []
for prefix in ["NC_PerGene_", "NC_GOBP_", "NC_HPO_", "CDS_GOBP_", "CDS_HPO_"]:
    df = 0
    ind = 1
    for file in file_list:
        if permute:
            v = pd.read_csv(out_dir + "/" + prefix + file.replace(prefix_run + ".ClosestVar.bed.gz", "ConvNearest.ForCombine.Prelim.Permute.Seed" + str(seed) + ".csv"))
        else:
            v = pd.read_csv(out_dir + "/" + prefix + file.replace(prefix_run + ".ClosestVar.bed.gz", "ConvNearest.ForCombine.Prelim.csv"))
        for index, row in v.iterrows():
            if row["Wilcoxon FDR"] < 0.25:
                out.append(list(row) + [prefix + "Conv", file.replace("." + prefix_run + ".Closest.bed.gz", "")])
        v = v[["Term", "Paired Mean Difference Corrected Quant", "Wilcoxon p-value Corrected Quant"]].set_index("Term")
        v.columns = [x + " " + file.replace("." + prefix_run + ".ClosestVar.bed.gz", "") for x in v.columns]
        if ind:
            df = v
            ind = 0
        else:
            df = df.join(v, how = "outer")
        
    cols_all = []
    for i in df.columns:
        if "Wilcoxon p-value" in i:
            cols_all.append(i)

    combined = []
    for index, row in df.iterrows():
        r = row[cols_all].dropna()
        #Change this when we get the rest of the data!
        if "PerGene" not in prefix:
            cut = len(file_list) - 1
        else:
            cut = len(file_list)//2
        if len(r) > cut:
            combined.append(combine_pvalues(r, method = "pearson")[1])
        else:
            combined.append(np.nan)
    df["Combined p-value all comps"] = combined
    df = df.dropna(subset = ["Combined p-value all comps"])
    df["Combined FDR all comps"] = fdrcorrection(df["Combined p-value all comps"])[1]
    df = df.sort_values("Combined p-value all comps")
    if permute:
        prefix = prefix + "Permute.Seed" + str(seed) + "_"
    df.to_csv(out_dir + "/Combined_Files/" + prefix + folder_name + "_Combined_ConvNearest.csv", sep = ",", index = True)
df = pd.DataFrame(out)
if permute:
    v = pd.read_csv(out_dir + "/" + prefix.replace("Permute.Seed" + str(seed) + "_", "") + file.replace(prefix_run + ".ClosestVar.bed.gz", "ConvNearest.ForCombine.Prelim.csv"))
else:
    v = pd.read_csv(out_dir + "/" + prefix + file.replace(prefix_run + ".ClosestVar.bed.gz", "ConvNearest.ForCombine.Prelim.csv"))
df.columns = list(v.columns) + ["VarType_Ontology", "Species Comparison"]
if permute:
    df.to_csv(out_dir + "/Combined_Files/" + folder_name + "_ConvNearest_AllFDRLessThan0.25_Permute.Seed" + str(seed) + ".csv", index = False)
else:
    df.to_csv(out_dir + "/Combined_Files/" + folder_name + "_ConvNearest_AllFDRLessThan0.25.csv", index = False)
    

#Combine p-values for the divergent sites
out = []
for prefix in ["NC_PerGene_", "NC_GOBP_", "NC_HPO_", "CDS_GOBP_", "CDS_HPO_"]:
    df = 0
    ind = 1
    for file in file_list:
        if permute:
            v = pd.read_csv(out_dir + "/" + prefix + file.replace(prefix_run + ".ClosestVar.bed.gz", "DivNearest.ForCombine.Prelim.Permute.Seed" + str(seed) + ".csv"))
        else:
            v = pd.read_csv(out_dir + "/" + prefix + file.replace(prefix_run + ".ClosestVar.bed.gz", "DivNearest.ForCombine.Prelim.csv"))
        for index, row in v.iterrows():
            if row["Wilcoxon FDR"] < 0.25:
                out.append(list(row) + [prefix + "Conv", file.replace("." + prefix_run + ".Closest.bed", "")])
        v = v[["Term", "Paired Mean Difference Corrected Quant", "Wilcoxon p-value Corrected Quant"]].set_index("Term")
        v.columns = [x + " " + file.replace("." + prefix_run + ".ClosestVar.bed.gz", "") for x in v.columns]
        if ind:
            df = v
            ind = 0
        else:
            df = df.join(v, how = "outer")
        
    cols_all = []
    for i in df.columns:
        if "Wilcoxon p-value" in i:
            cols_all.append(i)

    combined = []
    for index, row in df.iterrows():
        r = row[cols_all].dropna()
        #Change this when we get the rest of the data!
        if "PerGene" not in prefix:
            cut = len(file_list) - 1
        else:
            cut = len(file_list)//2
        if len(r) > cut:
            combined.append(combine_pvalues(r, method = "pearson")[1])
        else:
            combined.append(np.nan)
    df["Combined p-value all comps"] = combined
    df = df.dropna(subset = ["Combined p-value all comps"])
    df["Combined FDR all comps"] = fdrcorrection(df["Combined p-value all comps"])[1]
    df = df.sort_values("Combined p-value all comps")
    if permute:
        prefix = prefix + "Permute.Seed" + str(seed) + "_"
    df.to_csv(out_dir + "/Combined_Files/" + prefix + folder_name + "_Combined_DivNearest.csv", sep = ",", index = True)
df = pd.DataFrame(out)
if permute:
    v = pd.read_csv(out_dir + "/" + prefix.replace("Permute.Seed" + str(seed) + "_", "") + file.replace(prefix_run + ".ClosestVar.bed.gz", "DivNearest.ForCombine.Prelim.csv"))
else:
    v = pd.read_csv(out_dir + "/" + prefix + file.replace(prefix_run + ".ClosestVar.bed.gz", "DivNearest.ForCombine.Prelim.csv"))
df.columns = list(v.columns) + ["VarType_Ontology", "Species Comparison"]
if permute:
    df.to_csv(out_dir + "/Combined_Files/" + folder_name + "_DivNearest_AllFDRLessThan0.25_Permute.Seed" + str(seed) + ".csv", index = False)
else:
    df.to_csv(out_dir + "/Combined_Files/" + folder_name + "_DivNearest_AllFDRLessThan0.25.csv", index = False)
    

#Combine p-values for convergent and divergent
d_div = {}
d_conv = {}

for file in os.listdir(out_dir + "/Combined_Files"):
    if "Combined_DivNearest" in file:
        if permute:
            if "Permute.Seed" + str(seed) + "" in file:
                d_div["_".join(file.split("_")[0:2])] = file
        else:
            if "Permute" not in file:
                d_div["_".join(file.split("_")[0:2])] = file
    elif "Combined_ConvNearest" in file:
        if permute:
            if "Permute.Seed" + str(seed) + "" in file:
                d_conv["_".join(file.split("_")[0:2])] = file
        else:
            if "Permute" not in file:
                d_conv["_".join(file.split("_")[0:2])] = file 
print(d_conv, d_div)
for key in d_conv.keys():
    div = pd.read_csv(out_dir + "/Combined_Files/" + d_div[key]).set_index("Term")
    div.columns = [x + " Div" for x in div.columns]
    conv = pd.read_csv(out_dir + "/Combined_Files/" + d_conv[key]).set_index("Term")
    conv.columns = [x + " Conv" for x in conv.columns]
    combined = conv.join(div, how = "outer")
    
    combined_ps = []
    keep_cols = []
    for col in combined.columns:
        if "Combined p-value all comps" in col:
            keep_cols.append(col)
    for index, row in combined.iterrows():
        combined_ps.append(combine_pvalues(row[keep_cols].dropna(), method = "pearson")[1])
    combined["ConvDiv Combined p-value"] = combined_ps
    print("Convergent vs Divergent correlation: ", key, pearsonr(-np.log10(combined.dropna()[keep_cols[0]]), -np.log10(combined.dropna()[keep_cols[1]])), spearmanr(-np.log10(combined.dropna()[keep_cols[0]]), -np.log10(combined.dropna()[keep_cols[1]])))
    combined = combined.sort_values("ConvDiv Combined p-value")
    combined["ConvDiv Combined FDR"] = fdrcorrection(combined["ConvDiv Combined p-value"])[1]
    if permute:
        combined.to_csv(out_dir + "/Combined_Files/" + key + "_" + folder_name + "_ConvDiv_Combined_Permute.Seed" + str(seed) + ".csv")
    else:
        combined.to_csv(out_dir + "/Combined_Files/" + key + "_" + folder_name + "_ConvDiv_Combined.csv")