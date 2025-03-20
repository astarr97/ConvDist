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
from scipy.stats import norm
import os
import sys
from scipy.stats import combine_pvalues
from scipy.stats import binom_test

#hfont = {'fontname':'Arial'}
#plt.rcParams["font.family"] = "Arial"

#Code borrowed heavily from here: https://stackoverflow.com/questions/62375034/find-non-overlapping-area-between-two-kde-plots
#plt.rcParams.update(
#    {"text.usetex": False}
#)

file_list = sys.argv[1]
file_list = file_list.split(",")

folder_name = sys.argv[2]
spec_sup = float(sys.argv[3])
rem_single = bool(int(sys.argv[4]))
ref_species = sys.argv[5]
out_dir = folder_name + "_FWC_PWC_Results"

#Function to prepare to compute alpha
def prepare_alpha(fixed, poly, stat = "PhyloP447"):
    x = np.repeat("Fixed", len(list(fixed.index)))
    y = np.repeat("Polymorphic", len(list(poly.index)))
    poly["FixedOrPoly"] = y
    fixed["FixedOrPoly"] = x
    v2 = fixed[["FixedOrPoly", stat]]
    vv2 = poly[["FixedOrPoly", stat]]
    x2 = v2[stat]
    yvals2 = vv2[stat]
    vvv = pd.concat([v2, vv2])
    vvv["PhyloP"] = vvv[stat].astype(float)
    return vvv

#Code to compute alpha
def compute_alpha_new(vvv, dn_cut = 0.1, weighted = False, plot = False, sub = False, stat = "PhyloP", window = [-5, 12], title = ""):
    x0 = vvv[vvv["FixedOrPoly"].isin(["Fixed"])][stat].astype(float)
    x1 = vvv[vvv["FixedOrPoly"].isin(["Polymorphic"])][stat].astype(float)

    kde0 = gaussian_kde(x0, bw_method=0.3)
    kde1 = gaussian_kde(x1, bw_method=0.3)

    xmin = min(x0.min(), x1.min())
    xmax = max(x0.max(), x1.max())
    dx = 0.2 * (xmax - xmin) # add a 20% margin, as the kde is wider than the data
    xmin -= dx
    xmax += dx

    x = np.linspace(xmin, xmax, 500)
    kde0_x = kde0(x)
    kde1_x = kde1(x)

    to_find_inter = kde0_x - kde1_x
    prev = ""
    crosses = []
    cross_signs = []
    #If search start is 0, you get erroneous results for some pathological distributions
    #e.g., NonCoding, NRXN3 SpecSup250, MAF0.25, PhastCons > 0 
    #If search start is the median of the polymorphic distribution, then it is less stable
    search_start = np.median(x1)
    #search_start = 0
    for i in range(len(x)):
        if x[i] > search_start:
            if type(prev) == str:
                prev = to_find_inter[i]
            else:
                if prev*to_find_inter[i] <= 0:
                    crosses.append(i)
                    if to_find_inter[i] > 0:
                        cross_signs.append(1)
                    else:
                        cross_signs.append(0)
                prev = to_find_inter[i]
    crosses.append(499)
    alpha = 0
    
    start_compute = 0
    for i in range(len(crosses)):
        if i + 1 != len(crosses):
            if cross_signs[i]:
                start_compute = crosses[i]
                break
    if start_compute:
        #dn = np.trapz(kde0_x[start_compute:500], x[start_compute:500])
        #ds = np.trapz(kde0_x[0:start_compute], x[0:start_compute])
        #pn = np.trapz(kde1_x[start_compute:500], x[start_compute:500])
        #ps = np.trapz(kde1_x[0:start_compute], x[0:start_compute])
        
        dn = vvv[(vvv[stat] > x[start_compute]) & (vvv["FixedOrPoly"] == "Fixed")].shape[0]
        pn = vvv[(vvv[stat] > x[start_compute]) & (vvv["FixedOrPoly"] == "Polymorphic")].shape[0]
        ds = vvv[(vvv[stat] <= x[start_compute]) & (vvv["FixedOrPoly"] == "Fixed")].shape[0]
        ps = vvv[(vvv[stat] <= x[start_compute]) & (vvv["FixedOrPoly"] == "Polymorphic")].shape[0]
        
        
        
        #If this is changed, you accept less stable estimates leading to bad things when e.g. determining a confidence interval
        #However, when doing fisher's exact test I think it makes sense to change this to 0.01
        if dn != 0 and dn/(dn + ds) > dn_cut:
            alpha = 1 - (ds/dn)*(pn/ps)
        else:
            for i in range(len(x)):
                if x[i] <= search_start and x[i+1] >= search_start:
                    #print(i)
                    start_compute = i
            #dn = np.trapz(kde0_x[start_compute:500], x[start_compute:500])
            #ds = np.trapz(kde0_x[0:start_compute], x[0:start_compute])
            #pn = np.trapz(kde1_x[start_compute:500], x[start_compute:500])
            #ps = np.trapz(kde1_x[0:start_compute], x[0:start_compute])
            
            dn = vvv[(vvv[stat] > x[start_compute]) & (vvv["FixedOrPoly"] == "Fixed")].shape[0]
            pn = vvv[(vvv[stat] > x[start_compute]) & (vvv["FixedOrPoly"] == "Polymorphic")].shape[0]
            ds = vvv[(vvv[stat] <= x[start_compute]) & (vvv["FixedOrPoly"] == "Fixed")].shape[0]
            ps = vvv[(vvv[stat] <= x[start_compute]) & (vvv["FixedOrPoly"] == "Polymorphic")].shape[0]
            alpha = 1 - (ds/dn)*(pn/ps)
    else:
        for i in range(len(x)):
            if x[i] <= search_start and x[i+1] >= search_start:
                #print(i)
                start_compute = i

        #dn = np.trapz(kde0_x[start_compute:500], x[start_compute:500])
        #ds = np.trapz(kde0_x[0:start_compute], x[0:start_compute])
        #pn = np.trapz(kde1_x[start_compute:500], x[start_compute:500])
        #ps = np.trapz(kde1_x[0:start_compute], x[0:start_compute])
        
        dn = vvv[(vvv[stat] > x[start_compute]) & (vvv["FixedOrPoly"] == "Fixed")].shape[0]
        pn = vvv[(vvv[stat] > x[start_compute]) & (vvv["FixedOrPoly"] == "Polymorphic")].shape[0]
        ds = vvv[(vvv[stat] <= x[start_compute]) & (vvv["FixedOrPoly"] == "Fixed")].shape[0]
        ps = vvv[(vvv[stat] <= x[start_compute]) & (vvv["FixedOrPoly"] == "Polymorphic")].shape[0]
        alpha = 1 - (ds/dn)*(pn/ps)
        
    """if plot:
        fig, ax = plt.subplots(figsize=(10,6))
        ax.plot(x, kde0_x, color='b', label='Fixed')
        ax.fill_between(x, kde0_x, 0, color='b', alpha=0.2)
        ax.plot(x, kde1_x, color='orange', label='Polymorphic')
        ax.fill_between(x, kde1_x, 0, color='orange', alpha=0.2)
        handles, labels = plt.gca().get_legend_handles_labels()
        ax.legend(handles, labels, title='', frameon=False, fontsize = 12)
        #plt.xlabel("PhyloP", **hfont, fontsize = 18)
        ax.set_ylabel("Density", **hfont, fontsize = 18)
        ax.set_xlim(window[0], window[1])
        plt.tight_layout()
        if title:
            ax.set_title(title)
        #plt.show()"""
    tag = "dc > 0.1, likely stable"
    
    
    if dn == 0 or dn/(dn + ds) < 0.1:
        tag = "dc < 0.1, may be unstable"
    if plot:
        return alpha, x[start_compute], [x[i] for i in crosses], tag, [[dn, ds], [pn, ps]], fig
    else:
        return alpha, x[start_compute], [x[i] for i in crosses], tag, [[dn, ds], [pn, ps]]


def quantile_normalize(v_to_norm, vv_to_norm):
    #Requires that v always has fewer (or equal) number of rows
    assert(v_to_norm.shape[0] <= vv_to_norm.shape[0])
    
    #Upsample the lower one, maximum of one duplicate per entry added
    if (vv_to_norm.shape[0] - v_to_norm.shape[0])/v_to_norm.shape[0] > 1:
        v_samp = v_to_norm.sample(frac = (vv_to_norm.shape[0] - v_to_norm.shape[0])/v_to_norm.shape[0], replace = True)
    else:
        v_samp = v_to_norm.sample(frac = (vv_to_norm.shape[0] - v_to_norm.shape[0])/v_to_norm.shape[0], replace = False)
    
    #This index isn't used, just needed to allow the dataframes to concatenate
    v_samp.index = range(v_to_norm.shape[0], v_to_norm.shape[0] + v_samp.shape[0])
    v_to_norm_new = pd.concat([v_to_norm, v_samp]).sort_values(6)

    #Make new dataframes and join
    vv_to_norm_new = vv_to_norm.copy().sort_values(6)
    v_to_norm_new.columns = ["v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11", "v12"]
    vv_to_norm_new.columns = ["vv0", "vv1", "vv2", "vv3", "vv4", "vv5", "vv6", "vv7", "vv8", "vv9", "vv10", "vv11", "vv12"]
    
    #Re-index to join the dataframes
    v_to_norm_new.index = list(range(v_to_norm_new.shape[0]))
    vv_to_norm_new.index = list(range(vv_to_norm_new.shape[0]))
    v_quant = v_to_norm_new.join(vv_to_norm_new)

    #Make dataframe to normalize
    df = v_quant[["v6", "vv6"]].copy()
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    df.columns = ["v6 norm", "vv6 norm"]
    new_v1 = np.mean(df, axis = 1)
    new_v2 = np.mean(df, axis = 1)
    df["v6 norm"] = new_v1
    df["vv6 norm"] = new_v2
    return df[["v6 norm"]].join(v_to_norm_new), df[["vv6 norm"]].join(vv_to_norm_new)

def remove_singeltons(vv, folder, prefix, remove_multi=True, rel=False, ref_spec=ref_species):
    if rel:
        multi_foc = prefix + ".MultiSingelton.Focal.Rel.bed.gz"
        single_foc = prefix + ".Singelton.Focal.Rel.bed.gz"
    else:
        multi_foc = prefix + ".MultiSingelton.Focal.bed.gz"
        single_foc = prefix + ".Singelton.Focal.bed.gz"
    if ref_spec not in multi_foc:
        multi_foc = multi_foc.replace(".bed.gz", "." + ref_spec + ".bed.gz")
    if ref_spec not in single_foc:
        single_foc = single_foc.replace(".bed.gz", "." + ref_spec + ".bed.gz")
   
    if remove_multi:
        files = [single_foc, multi_foc]
    else:
        files = [single_foc]
    
    toss = []
    for file in files:
        x = pd.read_csv(folder + "/" + file, header = None, sep = "\t")
        toss = toss + [i[0] + ":" + str(i[1]) for i in zip(x[0], x[2])]
    vv["Chrom:Pos2"] = vv[0] + ":" + vv[2].astype(str)
    num_start = vv.shape[0]
    vv = vv[~vv["Chrom:Pos2"].isin(toss)].drop("Chrom:Pos2", axis = 1).copy()
    num_end = vv.shape[0]
    print("Proportion singeltons (all removed): ", 1-num_end/num_start)
    return vv
    
def read_fp(fixed_file, poly_file, spec_sup = 100, rem_single=rem_single):
    np.random.seed(6)
    v = pd.read_csv(fixed_file, sep = "\t", header = None)
    v = v[v[7] != "."]
    v = v[v[6] != "."]
    v[6] = v[6].astype(float)
    v = v[v[9] != "."]
    v[9] = v[9].astype(int)
    v = v[v[9] > spec_sup]
    
    vv = pd.read_csv(poly_file, sep = "\t", header = None)
    vv = vv[vv[7] != "."]
    vv = vv[vv[6] != "."]
    vv[6] = vv[6].astype(float)
    vv = vv[vv[9] != "."]
    vv[9] = vv[9].astype(int)
    vv = vv[vv[9] > spec_sup]
    
    v = v.sort_values(7).drop_duplicates([0, 2])
    vv = vv.sort_values(7).drop_duplicates([0, 2])
    
    vv = vv.sort_values(7).drop_duplicates([0, 2])
    
    if rem_single:
        folder = "../" + "_".join(fixed_file.split(".")[0].split("_")[0:2])
        prefix = fixed_file.split(".")[0]
        vv = remove_singeltons(vv, folder, prefix)

    v_nc = v[v[8] != 0].copy()
    vv_nc = vv[vv[8] != 0].copy()
    
    v_cds = v[v[8] == 0].copy()
    vv_cds = vv[vv[8] == 0].copy()
    
    #Reduce memory usage
    v = 0
    vv = 0
    
    v_cds["Corrected"] = v_cds[6] - (np.median(v_cds[6]) - np.median(vv_cds[6]))
    vv_cds["Corrected"] = vv_cds[6]

    v_cds["GEQ 0"] = np.maximum(v_cds[6], 0)
    vv_cds["GEQ 0"] = np.maximum(vv_cds[6], 0)

    v_cds["Corrected >= 0"] = v_cds["GEQ 0"] - (np.median(v_cds["GEQ 0"]) - np.median(vv_cds["GEQ 0"]))
    vv_cds["Corrected >= 0"] = vv_cds["GEQ 0"]
    
    v_nc["Corrected"] = v_nc[6] - (np.median(v_nc[6]) - np.median(vv_nc[6]))
    vv_nc["Corrected"] = vv_nc[6]

    v_nc["GEQ 0"] = np.maximum(v_nc[6], 0)
    vv_nc["GEQ 0"] = np.maximum(vv_nc[6], 0)

    v_nc["Corrected >= 0"] = v_nc["GEQ 0"] - (np.median(v_nc["GEQ 0"]) - np.median(vv_nc["GEQ 0"]))
    vv_nc["Corrected >= 0"] = vv_nc["GEQ 0"]
    
    print(v_cds)
    
    if len(v_cds.index) > len(vv_cds.index):
        vv_cds_new, v_cds_new = quantile_normalize(vv_cds, v_cds)

        v_cds_new = v_cds_new.sample(frac = 1, replace = False)
        v_cds_new = v_cds_new.drop_duplicates(["vv0", "vv2"])
        v_cds = v_cds_new.copy()
        v_cds = v_cds[["vv0", "vv2", "vv6", "vv7", "vv8", "vv9", "vv6 norm"]]
        v_cds.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        v_cds.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    
        vv_cds_new = vv_cds_new.sample(frac = 1, replace = False)
        vv_cds_new = vv_cds_new.drop_duplicates(["v0", "v2"])
        vv_cds = vv_cds_new.copy()
        vv_cds = vv_cds[["v0", "v2", "v6", "v7", "v8", "v9", "v6 norm"]]
        vv_cds.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        vv_cds.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    else:
        v_cds_new, vv_cds_new = quantile_normalize(v_cds, vv_cds)

        v_cds_new = v_cds_new.sample(frac = 1, replace = False)
        v_cds_new = v_cds_new.drop_duplicates(["v0", "v2"])
        v_cds = v_cds_new.copy()
        v_cds = v_cds[["v0", "v2", "v6", "v7", "v8", "v9", "v6 norm"]]
        v_cds.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        v_cds.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    
        vv_cds_new = vv_cds_new.sample(frac = 1, replace = False)
        vv_cds_new = vv_cds_new.drop_duplicates(["vv0", "vv2"])
        vv_cds = vv_cds_new.copy()
        vv_cds = vv_cds[["vv0", "vv2", "vv6", "vv7", "vv8", "vv9", "vv6 norm"]]
        vv_cds.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        vv_cds.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    
    if len(v_nc.index) > len(vv_nc.index):
        vv_nc_new, v_nc_new = quantile_normalize(vv_nc, v_nc)
        
        v_nc_new = v_nc_new.sample(frac = 1, replace = False)
        v_nc_new = v_nc_new.drop_duplicates(["vv0", "vv2"])
        v_nc = v_nc_new.copy()
        v_nc = v_nc[["vv0", "vv2", "vv6", "vv7", "vv8", "vv9", "vv6 norm"]]
        v_nc.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        v_nc.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    
        vv_nc_new = vv_nc_new.sample(frac = 1, replace = False)
        vv_nc_new = vv_nc_new.drop_duplicates(["v0", "v2"])
        vv_nc = vv_nc_new.copy()
        vv_nc = vv_nc[["v0", "v2", "v6", "v7", "v8", "v9", "v6 norm"]]
        vv_nc.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        vv_nc.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    else:
        v_nc_new, vv_nc_new = quantile_normalize(v_nc, vv_nc)
        v_nc_new = v_nc_new.sample(frac = 1, replace = False)
        v_nc_new = v_nc_new.drop_duplicates(["v0", "v2"])
        v_nc = v_nc_new.copy()
        v_nc = v_nc[["v0", "v2", "v6", "v7", "v8", "v9", "v6 norm"]]
        v_nc.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        v_nc.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    
        vv_nc_new = vv_nc_new.sample(frac = 1, replace = False)
        vv_nc_new = vv_nc_new.drop_duplicates(["vv0", "vv2"])
        vv_nc = vv_nc_new.copy()
        vv_nc = vv_nc[["vv0", "vv2", "vv6", "vv7", "vv8", "vv9", "vv6 norm"]]
        vv_nc.columns = [0, 2, 6, 7, 8, 9, "Corrected Quant"]
        vv_nc.columns = ["Chrom", "Position2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport", "Corrected Quant"]
    print(v_cds)
    print(vv_nc)
    return v_cds, vv_cds, v_nc, vv_nc


#Test whether this works for FWC PWC!
#Define function to search for convergent evolution at gene set level
def shuffle_fp(vvv):
    shuffle = vvv.copy()
    shuffled = list(shuffle["FixedOrPoly"].sample(frac = 1, replace = False))
    shuffle["FixedOrPoly"] = shuffled
    return shuffle

def run_gene_set_test(v, vv, gene_set, min_var_cut=25, stat_test="Corrected Quant", stat_print="PhyloP", permute = False):
    out = []
    out_g0 = []
    if type(gene_set) is dict:
        for key in gene_set.keys():
            v2 = v[v["NearestGene"].isin(gene_set[key])].copy()
            vv2 = vv[vv["NearestGene"].isin(gene_set[key])].copy()
            if len(v2.index) > min_var_cut and len(vv2.index) > min_var_cut:
                vvv = prepare_alpha(v2, vv2, stat = stat_test)
                alpha = compute_alpha_new(vvv, plot = False, dn_cut = 0.1)
                alpha_all = alpha[0]*alpha[4][0][0]/(alpha[4][0][0] + alpha[4][0][1])
                fp = fisher_exact(alpha[4], alternative = "greater")
                if permute:
                    shuffle_alphas = []
                    shuffle_alphas_all = []
                    for i in range(1000):
                        if i % 100 == 0:
                            print(i)
                        shuffled = shuffle_fp(vvv)
                        alpha_shuf = compute_alpha_new(shuffled, plot = False, dn_cut = 0.1)
                        shuffle_alphas.append(alpha_shuf[0])
                        shuffle_alphas_all.append(alpha_shuf[0]*alpha_shuf[4][0][0]/(alpha_shuf[4][0][0] + alpha_shuf[4][0][1]))
                    z = (alpha[0] - np.mean(shuffle_alphas))/np.std(shuffle_alphas)
                    z_all = (alpha_all - np.mean(shuffle_alphas_all))/np.std(shuffle_alphas_all)
                    permute_p = norm.sf(z)
                    permute_p_all = norm.sf(z_all)
                if permute:
                    out.append([key, np.median(v2[stat_print]), v2.shape[0], np.median(vv2[stat_print]), vv2.shape[0], mwu(v2[stat_test], vv2[stat_test], alternative = "greater")[1], fp[1], alpha[0], alpha[1], alpha[2], alpha[3], alpha[4], z, permute_p, np.median(shuffle_alphas), alpha_all, z_all, permute_p_all, np.median(shuffle_alphas_all)])
                else:
                    out.append([key, np.median(v2[stat_print]), v2.shape[0], np.median(vv2[stat_print]), vv2.shape[0], mwu(v2[stat_test], vv2[stat_test], alternative = "greater")[1], fp[1], alpha[0], alpha[1], alpha[2], alpha[3], alpha[4]])
    elif type(gene_set) is list:
        for gene in gene_set:
            v2 = v[v["NearestGene"].isin([gene])].copy()
            vv2 = vv[vv["NearestGene"].isin([gene])].copy()
            if len(v2.index) > min_var_cut and len(vv2.index) > min_var_cut:
                vvv = prepare_alpha(v2, vv2, stat = stat_test)
                try:
                    alpha = compute_alpha_new(vvv, plot = False, dn_cut = 0.1)
                    alpha_all = alpha[0]*alpha[4][0][0]/(alpha[4][0][0] + alpha[4][0][1])
                    fp = fisher_exact(alpha[4], alternative = "greater")
                    if permute:
                        shuffle_alphas = []
                        shuffle_alphas_all = []
                        for i in range(1000):
                            shuffled = shuffle_fp(vvv)
                            alpha_shuf = compute_alpha_new(shuffled, plot = False, dn_cut = 0.1)
                            shuffle_alphas.append(alpha_shuf[0])
                            shuffle_alphas_all.append(alpha_shuf[0]*alpha_shuf[4][0][0]/(alpha_shuf[4][0][0] + alpha_shuf[4][0][1]))
                        z = (alpha[0] - np.mean(shuffle_alphas))/np.std(shuffle_alphas)
                        z_all = (alpha_all - np.mean(shuffle_alphas_all))/np.std(shuffle_alphas_all)
                        permute_p = norm.sf(z)
                        permute_p_all = norm.sf(z_all)
                    if permute:
                        out.append([gene, np.median(v2[stat_print]), v2.shape[0], np.median(vv2[stat_print]), vv2.shape[0], mwu(v2[stat_test], vv2[stat_test], alternative = "greater")[1], fp[1], alpha[0], alpha[1], alpha[2], alpha[3], alpha[4], z, permute_p, np.median(shuffle_alphas), alpha_all, z_all, permute_p_all, np.median(shuffle_alphas_all)])
                    else:
                        out.append([gene, np.median(v2[stat_print]), v2.shape[0], np.median(vv2[stat_print]), vv2.shape[0], mwu(v2[stat_test], vv2[stat_test], alternative = "greater")[1], fp[1], alpha[0], alpha[1], alpha[2], alpha[3], alpha[4]])
                except:
                    print(gene)


    df = pd.DataFrame(out)
    if df.shape[0]:
        if permute:
            df["FDR"] = fdrcorrection(df[13])[1]
            df["FDR all"] = fdrcorrection(df[17])[1]
            df.columns = ["Term", "Median " + stat_print + " FWC", "Number FWC Sites", "Median " + stat_print + " PWC", "Number PWC Sites", "MWU p-value " + stat_test, "Fisher exact p-value " + stat_test, "Alpha", "Cutoff used", "Cutoffs", "Tag", "Table", "Permuted z-score", "Permuted p-value", "Median shuffled alpha",  "Alpha all", "Permuted z-score all", "Permuted p-value all", "Median shuffled alpha all", "Permuted FDR", "Permuted FDR all"]
        else:
            df["FDR"] = fdrcorrection(df[5])[1]
            df.columns = ["Term", "Median " + stat_print + " FWC", "Number FWC Sites", "Median " + stat_print + " PWC", "Number PWC Sites", "MWU p-value " + stat_test, "Fisher exact p-value " + stat_test, "Alpha", "Cutoff used", "Cutoffs", "Tag", "Table", "MWU FDR"]
    prop_fwc_pt_norm = []
    binom_pvalues = []
    prop_fwc_gw = len(v.index)/(len(v.index) + len(vv.index))
    for index, row in df.iterrows():
        prop_fwc_pt = row["Number FWC Sites"]/(row["Number FWC Sites"] + row["Number PWC Sites"])
        binom_p = binom_test(row["Number FWC Sites"], row["Number FWC Sites"] + row["Number PWC Sites"], p = prop_fwc_gw)
        prop_fwc_pt_norm.append(prop_fwc_pt/prop_fwc_gw)
        binom_pvalues.append(binom_p)
    df["Prop FWC PerTerm/Prop FWC Genome-Wide"] = prop_fwc_pt_norm
    df["Binomial p-value FWC vs PWC"] = binom_pvalues
    return df
    

gobp = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Ontologies/GOBP_AccelEvol_Input.txt", sep= "\t")
d_BP = {}

for index, row in gobp.iterrows():
    d_BP[row["Term"]] = row["Genes"].split(";")
    
hpo = pd.read_csv("/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Ontologies/HPO_AccelEvol_Input.txt", sep= "\t")
d_HPO = {}

for index, row in hpo.iterrows():
    d_HPO[row["Term"]] = row["Genes"].split(";")
    

if out_dir not in os.listdir():
    os.mkdir(out_dir)

for fixed_file in file_list:
    poly_file = fixed_file.replace("FiltConv", "FiltPoly")
    v_nc = 0
    vv_nc = 0
    v_cds = 0
    vv_cds = 0
    v_cds, vv_cds, v_nc, vv_nc = read_fp(fixed_file, poly_file)
    
    genes_cds = list(np.intersect1d(v_cds["NearestGene"], vv_cds["NearestGene"]))
    genes_nc = list(np.intersect1d(v_nc["NearestGene"], vv_nc["NearestGene"]))
    
    #Run the tests on HPO, GOBP, and per gene
    out_suffix = fixed_file.split(".")[0] + ".FWC_PWC.Prelim.csv"

    df = run_gene_set_test(v_nc, vv_nc, d_BP, min_var_cut = 50)
    df.sort_values("MWU p-value Corrected Quant").to_csv(out_dir + "/NC_GOBP_" + out_suffix, index = False)

    df = run_gene_set_test(v_nc, vv_nc, d_HPO, min_var_cut = 50)
    df.sort_values("MWU p-value Corrected Quant").to_csv(out_dir + "/NC_HPO_" + out_suffix, index = False)
    
    df = run_gene_set_test(v_nc, vv_nc, genes_nc, min_var_cut = 50)
    df.sort_values("MWU p-value Corrected Quant").to_csv(out_dir + "/NC_PerGene_" + out_suffix, index = False)

    df = run_gene_set_test(v_cds, vv_cds, d_BP, min_var_cut = 25)
    if df.shape[0]:
        df.sort_values("MWU p-value Corrected Quant").to_csv(out_dir + "/CDS_GOBP_" + out_suffix, index = False)

    df = run_gene_set_test(v_cds, vv_cds, d_HPO, min_var_cut = 25)
    if df.shape[0]:
        df.sort_values("MWU p-value Corrected Quant").to_csv(out_dir + "/CDS_HPO_" + out_suffix, index = False)
    
    df = run_gene_set_test(v_cds, vv_cds, genes_cds, min_var_cut = 10)
    if df.shape[0]:
        df.sort_values("MWU p-value Corrected Quant").to_csv(out_dir + "/CDS_PerGene_" + out_suffix, index = False)


if "Combined_Files" not in os.listdir(out_dir):
    os.mkdir(out_dir + "/Combined_Files")
out = []
for prefix in ["NC_PerGene_", "NC_GOBP_", "NC_HPO_", "CDS_GOBP_", "CDS_HPO_"]:

    df = 0
    ind = 1
    for fixed_file in file_list:
        v = pd.read_csv(out_dir + "/" + prefix + fixed_file.split(".")[0] + ".FWC_PWC.Prelim.csv")
        for index, row in v.iterrows():
            if row["MWU FDR"] < 0.25:
                out.append(list(row) + [prefix[:-1], fixed_file.split(".")[0]])
        v = v[["Term", "Alpha", "Fisher exact p-value Corrected Quant", "MWU p-value Corrected Quant"]].set_index("Term")
        v.columns = [x + " " + fixed_file.split(".")[0] for x in v.columns]
        if ind:
            df = v
            ind = 0
        else:
            df = df.join(v, how = "outer")
        
    cols_all = []
    for i in df.columns:
        if "MWU p-value" in i:
            cols_all.append(i)

    combined = []
    for index, row in df.iterrows():
        r = row[cols_all].dropna()
        cut = 1
        if len(r) > cut:
            combined.append(combine_pvalues(r, method = "pearson")[1])
        else:
            combined.append(np.nan)
    df["Combined MWU p-value all comps"] = combined
    df = df.dropna(subset = ["Combined MWU p-value all comps"])
    df["Combined MWU FDR all comps"] = fdrcorrection(df["Combined MWU p-value all comps"])[1]
    df = df.sort_values("Combined MWU p-value all comps")
    df.to_csv(out_dir + "/Combined_Files/" + prefix + folder_name + "_Combined_FWC_PWC.csv", sep = ",", index = True)
df = pd.DataFrame(out)
v = pd.read_csv(out_dir + "/" + prefix + fixed_file.split(".")[0] + ".FWC_PWC.Prelim.csv")
df.columns = list(v.columns) + ["VarType_Ontology", "Species"]
df.to_csv(out_dir + "/Combined_Files/" + folder_name + "_FWC_PWC_AllFDRLessThan0.25.csv", index = False)