import sys
import os

file = sys.argv[1]
focal_species = sys.argv[2]
other_species = sys.argv[3]
out_dir = sys.argv[4]

com_snp = "/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halSnps --minSpeciesForSnp 1 --refSequence REFSEQ --start START --length 1 --tsv TSV_NAME /scratch/users/astarr97/PhyloP/hg38.447way.hal " + focal_species + " " + other_species + "\n"

o = open(file)
os.mkdir(out_dir + "/All")
os.mkdir(out_dir + "/run1")
out = open(out_dir + "/run1/" + "run1.sh", 'w')
out.write("#!/bin/bash\n#SBATCH --time=36:00:00\n#SBATCH -p hns,hbfraser\n#SBATCH --mem=8GB\n\n")

c = 0
tot = 1
for line in o:
    c += 1
    l = line.split("\t")
    
    out.write(com_snp.replace("REFSEQ", l[0]).replace("START", l[1]).replace("TSV_NAME", l[0] + "_" + l[1] + ".tsv"))
    if c % 10000 == 0:
        out.write("\n")
        out.write("cat *.tsv > Run" + str(tot) + ".tsv\n")
        #out.write("python ../../remHead_callConvDiv.py Run" + str(tot) + ".tsv " + ignore_first + "\n")
        out.write("mv Run" + str(tot) + ".tsv " + "../All\n")
        out.write("rm *.tsv\n")
        out.close()
        tot += 1
        os.mkdir(out_dir + "/run" + str(tot))
        out = open(out_dir + "/run" + str(tot) + "/run" + str(tot) + ".sh", 'w')
        out.write("#!/bin/bash\n#SBATCH --time=36:00:00\n#SBATCH -p hns,hbfraser\n#SBATCH --mem=8GB\n\n")

out.write("cat *.tsv > Run" + str(tot) + ".tsv\n")
#out.write("python ../../remHead_callConvDiv.py Run" + str(tot) + ".tsv " + ignore_first + "\n")
out.write("mv Run" + str(tot) + ".tsv " + "../All\n")
out.write("rm *.tsv\n")
out.close()
