import sys

file = sys.argv[1]
assert(file.startswith("ToMakeConfig_"))

o = open(file)
out = open(file.replace("ToMakeConfig_", "Config_GetDL_"), 'w')

k = []
for line in o:
    k.append(line.replace("\n", ""))
    
for j in range(len(k)):
    for i in range(len(k)):
        if j < i:
            out.write(k[j] + "," + k[i] + "\n")
o.close()
out.close()
