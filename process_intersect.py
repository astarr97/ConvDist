import sys

file = sys.argv[1]

assert(".Intersect.bed" in file)

o = open(file)
out = open(file.replace(".Intersect.bed", ".Ready.bed"), 'w')

c = 1
for line in o:
    l = line.split("\t")
    if l[10] != ".":
        out.write("\t".join(l[0:10] + l[13:16]) + "\n")
o.close()
out.close()
