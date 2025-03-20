import sys

file = sys.argv[1]

assert(".ToFix.bed" in file)

o = open(file)
out = open(file.replace(".ToFix.bed", ".bed"), 'w')

for line in o:
    l = line.split("\t")
    out.write("\t".join(l[0:6] + [l[9], l[14], l[15], l[19]]) + "\n")
o.close()
out.close()
