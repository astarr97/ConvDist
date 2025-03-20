import sys

file = sys.argv[1]

assert("Closest.bed" in file)

o = open(file)
out = open(file.replace("Closest.bed", "ReClose.bed"), 'w')

c = 1
for line in o:
    l = line.replace("\n", "").split("\t")
    if l[1] != l[26]:
        out.write("\t".join([l[0], l[1], l[2]] + l[6:10] + l[16:24] + l[25:]) + "\n")
o.close()
out.close()
