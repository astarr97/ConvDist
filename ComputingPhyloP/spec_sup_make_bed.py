import sys

file = sys.argv[1]

assert(".txt" in file)

o = open(file)
out = open(file.replace(".txt", ".bed"), 'w')

for line in o:
    l = line.split("\t")
    out.write("\t".join([l[0], str(int(l[1])), str(int(l[1]) + int(l[2])), l[4]]))
out.close()
o.close()
