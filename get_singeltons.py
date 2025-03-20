import sys

file = sys.argv[1]
out_prefix = sys.argv[2]

assert(".bed" in file)

o = open(file)
out_foc_s = open(out_prefix + ".Singelton.Focal.ToFix.bed", 'w')
out_rel_s = open(out_prefix + ".Singelton.Rel.ToFix.bed", 'w')
out_out_s = open(out_prefix + ".Singelton.Out.ToFix.bed", 'w')
out_foc_ms = open(out_prefix + ".MultiSingelton.Focal.ToFix.bed", 'w')
out_rel_ms = open(out_prefix + ".MultiSingelton.Rel.ToFix.bed", 'w')
out_out_ms = open(out_prefix + ".MultiSingelton.Out.ToFix.bed", 'w')

for line in o:
    l = line.replace("\n", "").split("\t")
    if l[16] == "S":
        out_foc_s.write("\t".join(l[0:3]) + "\n")
    if l[17] == "S":
        out_rel_s.write("\t".join(l[0:3]) + "\n")
    if l[18] == "S":
        out_out_s.write("\t".join(l[0:3]) + "\n")
    if l[16] == "MS":
        out_foc_ms.write("\t".join(l[0:3]) + "\n")
    if l[17] == "MS":
        out_rel_ms.write("\t".join(l[0:3]) + "\n")
    if l[18] == "MS":
        out_out_ms.write("\t".join(l[0:3]) + "\n")
o.close()
out_foc_s.close()
out_rel_s.close()
out_out_s.close()
out_foc_ms.close()
out_rel_ms.close()
out_out_ms.close()
