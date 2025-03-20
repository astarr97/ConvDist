import sys
import numpy as np
file = sys.argv[1]
ignore_first = int(sys.argv[2])

assert(".tsv" in file)

o = open(file)
out = open(file.replace(".tsv", ".bed"), 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    #The first species (after the chrom and position) is the reference species but the reference species is not in the comparison, so we remove it
    if ignore_first == 1:
        l = l[0:2] + l[3:]
    #The first species (after the chrom and position) is the reference species but we are doing rel so we want the first species to actually be the third species and move everything over
    elif ignore_first == 2:
        old4 = l[4]
        l[4] = l[2]
        l[2] = l[3]
        l[3] = old4
    if l[0] != "refSequence":
        if len(np.unique([x.upper() for x in l[4:]])) == 1:
            anc_same = "AncestralSame"
        else:
            anc_same = "AncestralNotSame"
        if l[2].upper() == l[3].upper():
            conv_or_div = "Convergent"
        else:
            conv_or_div = "Divergent"
        out.write("\t".join([l[0].split(".")[0], l[1], str(int(l[1]) + 1)] + [x.upper() for x in l[2:]] + [conv_or_div, anc_same]) + "\n")