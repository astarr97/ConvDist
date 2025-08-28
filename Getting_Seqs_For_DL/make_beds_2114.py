import sys
import pandas as pd

file = sys.argv[1]
species = sys.argv[2]
rel = int(sys.argv[3])
v = pd.read_csv(file, sep = "\t", header = None)

if not rel:
    #7 is the focal species no. 1 base, 8 is the focal species no. 2, 9 is the related species no. 1 base, 11 is the related species no. 2 base
    #We need the related species because there is not a guarantee that the ancestral sequence will match the related species due to the exclusion of some species from the initial variant pulling
    #For example, we use the LCA of bovidae, cetaceans, and hippos for studying cetaceans, but we do not do any filtering on hippo
    #As a result, the orca base can appear to be the ancestral base 
    v["Position1"] = v[0] + ":" + v[2].astype(str) + ";" + v[7] + ";" + v[8] + ";" + v[9] + ";" + v[11]
    
    #18 is the focal species no. 1 closest site base, 19 is the related species no. 1 closest site base
    v["Position2"] = v[15] + ":" + v[17].astype(str) + ";" + v[18] + ";" + v[19]
    
    #29 is the focal species no. 2 closest site base, 30 is the related species no. 2 closest site base
    v["Position3"] = v[26] + ":" + v[28].astype(str) + ";" + v[29] + ";" + v[30]
else:
    print("Should be Rel")
    #Even though we are doing rel, we don't need to swap the focal and related species for this due to how it is formatted
    v["Position1"] = v[0] + ":" + v[2].astype(str) + ";" + v[7] + ";" + v[8] + ";" + v[9] + ";" + v[11]
    
    #The Loxodonta_africana and Trichechus_manatus Rel files are specially formatted because of the weird uncertainty with hyraxes
    #As a result, this is unchanged from the non-Rel
    if "Loxodonta_africana" == file.split(".")[0] or "Trichechus_manatus" == file.split(".")[0]:
        v["Position2"] = v[15] + ":" + v[17].astype(str) + ";" + v[18] + ";" + v[19]
    else:
        v["Position2"] = v[15] + ":" + v[17].astype(str) + ";" + v[19] + ";" + v[18]
    #The Loxodonta_africana and Trichechus_manatus Rel files are specially formatted because of the weird uncertainty with hyraxes
    #As a result, this is unchanged from the non-Rel
    if "Loxodonta_africana" == file.split(".")[1] or "Trichechus_manatus" == file.split(".")[1]:
        v["Position3"] = v[26] + ":" + v[28].astype(str) + ";" + v[29] + ";" + v[30]
    else:
        v["Position3"] = v[26] + ":" + v[28].astype(str) + ";" + v[30] + ";" + v[29]

#This whole annoying code block is all about adding the version number back to the chromosomes...
#If orca, we need to add the ".1" to the end
if species == "Orcinus_orca":
    v[0] = v[0] + ".1"
    v[15] = v[15] + ".1"
    v[26] = v[26] + ".1"
    
### NEED TO TEST THAT THIS WORKS ###
elif species == "Homo_sapiens":
    new = []
    for index, l in v.iterrows():
        if l[0] in ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]:
            pass
        else:
            if l[0] in ["GL000008", "GL000009", "GL000205", "GL000216"]:
                l[0] = l[0] + ".2"
            else:
                l[0] = l[0] + ".1"
        if l[15] in ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]:
            pass
        else:
            if l[15] in ["GL000008", "GL000009", "GL000205", "GL000216"]:
                l[15] = l[15] + ".2"
            else:
                l[15] = l[15] + ".1"
        if l[26] in ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]:
            pass
        else:
            if l[26] in ["GL000008", "GL000009", "GL000205", "GL000216"]:
                l[26] = l[26] + ".2"
            else:
                l[26] = l[26] + ".1"
        new.append(l)
    v = pd.DataFrame(l)
elif species == "CanFam4":
    if l[0] in ["chr" + str(i) for i in range(1, 39)] + ["chrX", "chrY", "chrM"]:
        pass
    else:
        if not l[0].endswith("v1"):
            l[0] = l[0] + "v1"
        else:
            pass
    if l[15] in ["chr" + str(i) for i in range(1, 39)] + ["chrX", "chrY", "chrM"]:
        pass
    else:
        if not l[15].endswith("v1"):
            l[15] = l[15] + "v1"
        else:
            pass
    if l[26] in ["chr" + str(i) for i in range(1, 39)] + ["chrX", "chrY", "chrM"]:
        pass
    else:
        if not l[26].endswith("v1"):
            l[26] = l[26] + "v1"
        else:
            pass
elif species == "Mus_musculus":  
        point2 = ["GL456021", "GL456033", "GL455992", "GL455999", "GL456045", "GL456028", "JH792831", "GL456049", "GL456002", "GL456012", "GL456024", "GL456003", "KQ030493", "GL456054", "GL456053", "GL456022", "GL456017", "GL456001", "GL456060", "GL456042", "GL456026", "KB469741", "GL456044"]
        point3 = ["KB469738"]
        if l[0] in ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY", "chrM"] or l[0].endswith("_random") or l[0].startswith("chrUn"):
            pass
        elif l[0] in point2:
            l[0] = l[0] + ".2"
        elif l[0] in point3:
            l[0] = l[0] + ".3"
        else:
            l[0] = l[0] + ".1"
            
        if l[15] in ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY", "chrM"] or l[15].endswith("_random") or l[15].startswith("chrUn"):
            pass
        elif l[15] in point2:
            l[15] = l[15] + ".2"
        elif l[15] in point3:
            l[15] = l[15] + ".3"
        else:
            l[15] = l[15] + ".1"
            
        if l[26] in ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY", "chrM"] or l[26].endswith("_random") or l[26].startswith("chrUn"):
            pass
        elif l[26] in point2:
            l[26] = l[26] + ".2"
        elif l[26] in point3:
            l[26] = l[26] + ".3"
        else:
            l[26] = l[26] + ".1"
elif species == "Nycticebus_pygmaeus":
    pass
else:
    v[0] = v[0] + ".1"
    v[15] = v[15] + ".1"
    v[26] = v[26] + ".1"


#Write out the positions
v[[0, 1, 2, "Position1"]].dropna().to_csv(file.replace("ClosestVar.bed.gz", "Focal.bed"), sep = "\t", header = None, index = False)
v[[15, 16, 17, "Position2"]].dropna().to_csv(file.replace("ClosestVar.bed.gz", "Closest_S1.bed"), sep = "\t", header = None, index = False)
v[[26, 27, 28, "Position3"]].dropna().to_csv(file.replace("ClosestVar.bed.gz", "Closest_S2.bed"), sep = "\t", header = None, index = False)
