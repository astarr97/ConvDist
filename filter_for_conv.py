import sys

file = sys.argv[1]

assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".FiltConv.bed"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    #Get focal species-derived changes
    if l[6] == "Fo":
        #That are fixed within the clade
        if l[7] == "F":
            #That have at least 1/4 of the fixed species in the alignment
            if float(l[8]) > 0.25:
                #And that are not equal to any of the relative or outgroup species
                if l[13] == "0" and l[14] == "0":
                    out.write("\t".join(l[0:6]) + "\n")
o.close()
out.close()

o = open(file)
out = open(file.replace(".bed", ".FiltConv.Rel.bed"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    #Get related species-derived changes
    if l[6] == "Re":
        #That are fixed within the clade
        if l[9] == "F":
            #That have at least 1/4 of the fixed species in the alignment
            if float(l[10]) > 0.25:
                #And that are not equal to any of the focal or outgroup species
                if l[13] == "0" and l[15] == "0":
                    out.write("\t".join(l[0:6]) + "\n")
o.close()
out.close()
