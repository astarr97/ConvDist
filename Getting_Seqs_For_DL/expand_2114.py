import sys

#First argument is bed file, second argument is file containing contig lengths for the reference genome
file = sys.argv[1]
contigs = sys.argv[2]

assert(".bed" in file)

d = {}

#Read in the contig lengths and create a dictionary
o = open(contigs)

for line in o:
    l = line.replace("\n", "").split("\t")
    d[l[0]] = int(l[2])

o = open(file)
out = open(file.replace(".bed", ".2114.bed"), 'w')
for line in o:
    #Expand region to be 2114 bases
    l = line.split("\t")
    l[1] = str(int(l[1]) - 1056)
    l[2] = str(int(l[2]) + 1057)
    
    #If this is a valid region based on the size of the contig, write out
    if int(l[1]) >= 0 and int(l[2]) < d[l[0]]:
        out.write("\t".join(l))
o.close()
out.close()
