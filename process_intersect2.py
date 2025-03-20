import sys

file = sys.argv[1]
species_list = sys.argv[2].split(",")

assert("ToFix.bed" in file)

o = open(file)
out = open(file.replace("ToFix.bed", "ConvDiv.txt"), 'w')

out.write("\t".join(["Chrom", "Pos1", "Pos2", "PhyloP", "NearestGene", "NearestDist", "SpeciesSupport"] + species_list + ["ConvOrDiv", "AncSame"]) + "\n")
for line in o:
    l = line.replace("\n", "").split("\t")
    out.write("\t".join([l[0], l[1], l[2]] + l[6:10] + l[16:24]) + "\n")
    
o.close()
out.close()
