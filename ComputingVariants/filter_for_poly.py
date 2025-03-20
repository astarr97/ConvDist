import sys

file = sys.argv[1]

assert(".bed" in file)

o = open(file)
out = open(file.replace(".bed", ".FiltPoly.bed"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    #Get focal species-derived changes
    if l[6] == "Fo":
        #That are polymorphic within the clade
        if l[7] == "P":
            #That have at least 1/4 of the fixed species in the alignment
            if float(l[8]) > 0.25:
                #And that are equal to any of the relative or outgroup species
                if l[13] == "1" or l[14] == "1":
                    out.write("\t".join(l[0:6]) + "\n")
o.close()
out.close()

#Repeat for the related species
o = open(file)
out = open(file.replace(".bed", ".FiltPoly.Rel.bed"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    #Get related species-derived changes
    if l[6] == "Re":
        #That are polymorphic within the clade
        if l[9] == "P":
            #That have at least 1/4 of the fixed species in the alignment
            if float(l[10]) > 0.25:
                #And that are equal to any of the focal or outgroup species
                if l[13] == "1" or l[15] == "1":
                    out.write("\t".join(l[0:6]) + "\n")
o.close()
out.close()
