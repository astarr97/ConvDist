o = open("geneAnnotation_CDS_Pteropus_alecto.sort.bed")
out = open("geneAnnotation_CDS_Pteropus_alecto.sort.fixed.bed", 'w')

for line in o:
    l = line.split("\t")
    l[0] = l[0] + ".1"
    out.write("\t".join(l))
o.close()
out.close()
