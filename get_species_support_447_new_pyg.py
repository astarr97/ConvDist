import sys
from collections import Counter

in_maf = sys.argv[1]
out_file = sys.argv[2]
ref_species = sys.argv[3]

o = open(in_maf)
out = open(out_file, 'w')
spec_count = 0
information = []
for line in o:
    if "#" not in line:
        if line == "a\n" and spec_count:
            out.write("\t".join(information + [str(spec_count)]) + "\n")
            spec_count = 0
        else:
            
            if line[0] == "s":
                string = line.split()[6]
                counts = Counter(string)
                if "N" not in string and counts["-"]/len(string) < 0.1:
                    spec_count += 1
            if ref_species in line:
                l = line.replace("\n", "").split()
                information = [".".join(l[1].split(".")[1:]), l[2], l[3], l[4]]
out.write("\t".join(information + [str(spec_count)]) + "\n")
