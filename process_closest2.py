import sys

file = sys.argv[1]

assert("Closest.bed" in file)

o = open(file)
out = open(file.replace("Closest.bed", "ClosReady.bed").replace("Closest", ""), 'w')

c = 1
for line in o:
    l = line.replace("\n", "").split("\t")
    if l[1] != l[27]:
        out.write(line)
o.close()
out.close()
