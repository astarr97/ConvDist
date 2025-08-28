import pandas as pd
import sys

bed_file = sys.argv[1]
tsv_file = sys.argv[2]
spec = sys.argv[3]

assert(bed_file.endswith(".bed"))
assert(tsv_file.endswith(".tsv"))

#There should be 5 nodes in the list if it is not the closest sites, 3 if it is the closest sites
if "Focal" in bed_file:
    assert(len(spec.split(",")) == 5)
elif "Closest" in bed_file:
    assert(len(spec.split(",")) == 3)
else:
    print("Must have the string 'Focal' or the string 'Closest' in bed_file name")
    assert(False)

rev = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}

#Read in the bed file and tsv file, set up to join
v = pd.read_csv(bed_file, sep = "\t", header = None)
v.index = v[0] + ":" + v[1].astype(str)
tv = pd.read_csv(tsv_file, sep = "\t")
tv = tv[tv["refSequence"] != "refSequence"].copy()
tv.index = tv["refSequence"] + ":" + tv["refPosition"].astype(str)

#Remove duplicates and join

v = v.drop_duplicates([3], keep = False)
print(v.shape[0])
v = v.join(tv)
print(v.shape[0])
print(v.sort_values(spec.split(",")[1])[["refSequence", "refPosition", spec.split(",")[0], spec.split(",")[1]]])

if "Focal" in bed_file:
    #Get rid of sites that we lack enough information
    v = v.dropna(subset=["refSequence", "refPosition", spec.split(",")[0], spec.split(",")[1], spec.split(",")[2], spec.split(",")[3]])
    print(v.shape[0])
    out = []
    for index, row in v.iterrows():
        #Change the bases to the bases from the TSV for each species
        #print(row[3])
        l = row[3].split(";")
        l[1] = row[spec.split(",")[0]].upper()
        l[2] = row[spec.split(",")[1]].upper()
        l[3] = row[spec.split(",")[2]].upper()
        l[4] = row[spec.split(",")[3]].upper()
        row[3] = ";".join(l) + ";" + row[spec.split(",")[4]].upper()
        #print(row[3])
        out.append(row)
elif "Closest" in bed_file:
    v = v.dropna(subset=["refSequence", "refPosition", spec.split(",")[0], spec.split(",")[1]])
    print(v.shape[0])
    out = []
    for index, row in v.iterrows():
        #Change the bases to the bases from the TSV for each species
        #print(row[3])
        l = row[3].split(";")
        l[1] = row[spec.split(",")[0]].upper()
        l[2] = row[spec.split(",")[1]].upper()
        row[3] = ";".join(l) + ";" + row[spec.split(",")[2]].upper()
        #print(row[3])
        out.append(row)

#Create dataframe and write out
vn = pd.DataFrame(out)
vn = vn[[0, 1, 2, 3]].copy()
vn.to_csv(bed_file.replace(".bed", ".Comped.bed"), sep = "\t", header = False, index = False)
