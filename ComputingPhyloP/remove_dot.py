import sys

file = sys.argv[1]

assert(".ToFix.bed" in file)

o = open(file)
out = open(file.replace(".ToFix.bed", ".bed"), 'w')

for line in o:
    l = line.split("\t")
    l[0] = l[0].split(".")[0]
    out.write("\t".join(l))
o.close()
out.close()
