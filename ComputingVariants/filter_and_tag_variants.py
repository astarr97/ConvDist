import pandas as pd
import sys
import numpy as np
from collections import Counter

file = sys.argv[1]
spec_focus = sys.argv[2]
spec_rel = sys.argv[3]
spec_out = sys.argv[4]
if sys.argv[5] != "0":
    spec_focus_use = sys.argv[5].split(",")
else:
    spec_focus_use = []
if sys.argv[6] != "0":
    spec_rel_use = sys.argv[6].split(",")
else:
    spec_rel_use = []
if sys.argv[7] != "0":
    spec_out_use = sys.argv[7].split(",")
else:
    spec_out_use = []

assert(".tsv" in file)

#Determines if there is only one unique value in a 1-D series (to determine if fixed or polymorphic)
def is_unique(s):
    a = s.to_numpy().astype(str) # s.values (pandas<0.24)
    a = np.char.replace(a, "N", a[0])
    return (a[0] == a).all()

#Counts the proportion of valid bases in the file
def prop_base(s):
    counts = Counter()
    counts["N"] = 0
    counts["-"] = 0
    counts["A"] = 0
    counts["C"] = 0
    counts["G"] = 0
    counts["T"] = 0
    counts = Counter(s) + counts
    return np.round((counts["A"] + counts["T"] + counts["C"] + counts["G"])/len(s), 3)

#Determines if the intersection of two subsets (determined by group1 and group2) of the 1-D series s is non-null
#If it is non-null, then something in group1 equals something in group2
def nonnull_intersect(s, group1, group2):
    x = np.intersect1d(s[group1].replace("N", np.nan).dropna(), s[group2].replace("N", np.nan).dropna())
    if len(x):
        return 1
    else:
        return 0

#Quick note on the behavior of this function
#Returns NR if a site is fixed within the clade (e.g. counts are "A":5)
#Returns N if a site is not a singelton variant (e.g. counts are "A":3, "C":2)
#Returns MS if a site only has singelton variants but there are multiple (e.g. counts are "A":3, "C":1, "T":1)
#Returns S if a site only has one singelton variants (e.g. counts are "A":4, "C":1)
def is_singelton_poly(s):
    counts = Counter(s.replace("N", np.nan).dropna())
    if len(counts.keys()) == 1:
        return "NR"
    else:
        if len([i for i in counts.values() if i >= 2]) > 1:
            return "N"
        else:
            if len(counts.values()) > 2:
                return "MS"
            else:
                return "S"

#Decide who is derived based solely on the focal species and outgroup species
def simple_derived(s):
    if s[spec_rel] == s[spec_out]:
        return "Fo"
    elif s[spec_focus] == s[spec_out]:
        return "Re"
    else:
        return "Am"
    
#Read in the file
v = pd.read_csv(file, sep = "\t")
v = v.replace(np.nan, "N")
v = v[v["refSequence"] != "refSequence"]
v["refPosition"] = v["refPosition"].astype(int)

#Convert all variants to upper case
for col in list(v.columns):
    if "refSequence" != col and "refPosition" != col:
        v[col] = v[col].str.upper()
        
#Write the file out again, replacing it now that everything is uppercase
v.to_csv(file, sep = "\t", index = False)

#Filter out sites that are N or the same in spec_focus and spec_rel
v = v[(v[spec_focus] != v[spec_rel]) & (v[spec_focus].isin(["A", "C", "G", "T"])) & (v[spec_rel].isin(["A", "C", "G", "T"])) & (v[spec_out].isin(["A", "C", "G", "T"]))]

#First, add a flag for whether spec_foc is derived ("Fo"), spec_rel is derived ("Re"), or it is ambiguous "Am"
#This is solely based on the first three species
v["Simple_Derived"] = v.apply(simple_derived, axis = 1)

#Determine if a site is polymorphic within the focal species, "F" is fixed "P" is polymorphic
#Also write out the proportion of non-N sites
v["Focus_FixedOrPoly"] = v[[spec_focus] + spec_focus_use].apply(is_unique, axis = 1).replace(True, "F").replace(False, "P")
v["Focus_PropBase"] = v[[spec_focus] + spec_focus_use].apply(prop_base, axis = 1)

#Repeat this for the rel species and out species
v["Related_FixedOrPoly"] = v[[spec_rel] + spec_rel_use].apply(is_unique, axis = 1).replace(True, "F").replace(False, "P")
v["Related_PropBase"] = v[[spec_rel] + spec_rel_use].apply(prop_base, axis = 1)

v["Outgroup_FixedOrPoly"] = v[[spec_out] + spec_out_use].apply(is_unique, axis = 1).replace(True, "F").replace(False, "P")
v["Outgroup_PropBase"] = v[[spec_out] + spec_out_use].apply(prop_base, axis = 1)

#Also create flags for if any of the focal species have the same base as the rel species 
#And another if any of the focal species have the same base as the out species
v["FocusEqualRelated"] = v.apply(nonnull_intersect, args=([spec_focus] + spec_focus_use, [spec_rel] + spec_rel_use), axis = 1)
v["FocusEqualOutgroup"] = v.apply(nonnull_intersect, args=([spec_focus] + spec_focus_use, [spec_out] + spec_out_use), axis = 1)
v["RelatedEqualOutgroup"] = v.apply(nonnull_intersect, args=([spec_rel] + spec_rel_use, [spec_out] + spec_out_use), axis = 1)

#Finally, create a flag if the focal SNP is a singelton (S if it is, N if it is not, NR if it is not polymorphic so doesn't matter)
v["FocusIsSingelton"] = v[[spec_focus] + spec_focus_use].apply(is_singelton_poly, axis = 1)

#Create a flag if the relative SNP is a singelton
v["RelatedIsSingelton"] = v[[spec_rel] + spec_rel_use].apply(is_singelton_poly, axis = 1)

#Create a flag if the outgroup SNP is a singelton
v["OutgroupIsSingelton"] = v[[spec_out] + spec_out_use].apply(is_singelton_poly, axis = 1)


columns_keep = ["Simple_Derived", "Focus_FixedOrPoly", "Focus_PropBase", "Related_FixedOrPoly", "Related_PropBase", "Outgroup_FixedOrPoly", "Outgroup_PropBase", "FocusEqualRelated", "FocusEqualOutgroup", "RelatedEqualOutgroup", "FocusIsSingelton", "RelatedIsSingelton", "OutgroupIsSingelton"]

#Convert to bed format
v["refPosition1"] = v["refPosition"] + 1
v = v[["refSequence", "refPosition", "refPosition1", spec_focus, spec_rel, spec_out] + columns_keep]
v = v.dropna()

#Write out file ready for intersection with PhyloP
v.to_csv(file.replace(".tsv", "_ForPhyloP.bed"), sep = "\t", header = False, index = False)

