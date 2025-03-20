import sys
import os

file = sys.argv[1]
focal_species = sys.argv[2]
other_species = sys.argv[3]
out_dir = sys.argv[4]
ignore_first = sys.argv[5]

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
    #Hard coded to deal with species that don't all end in ".1"
    if focal_species == "Homo_sapiens":
        if l[0] in ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]:
            out.write(com_snp.replace("REFSEQ", l[0]).replace("START", l[1]).replace("TSV_NAME", l[0] + "_" + l[1] + ".tsv"))
        else:
            if l[0] in ["GL000008", "GL000009", "GL000205", "GL000216"]:
                out.write(com_snp.replace("REFSEQ", l[0] + ".2").replace("START", l[1]).replace("TSV_NAME", l[0] + ".2_" + l[1] + ".tsv"))
            else:
                out.write(com_snp.replace("REFSEQ", l[0] + ".1").replace("START", l[1]).replace("TSV_NAME", l[0] + ".1_" + l[1] + ".tsv"))
                
    #Version to deal with CanFam4
    elif focal_species == "CanFam4":
        if l[0] in ["chr" + str(i) for i in range(1, 39)] + ["chrX", "chrY", "chrM"]:
            out.write(com_snp.replace("REFSEQ", l[0]).replace("START", l[1]).replace("TSV_NAME", l[0] + "_" + l[1] + ".tsv"))
        else:
            if not l[0].endswith("v1"):
                out.write(com_snp.replace("REFSEQ", l[0] + "v1").replace("START", l[1]).replace("TSV_NAME", l[0] + "v1_" + l[1] + ".tsv"))
            else:
                out.write(com_snp.replace("REFSEQ", l[0]).replace("START", l[1]).replace("TSV_NAME", l[0] + "_" + l[1] + ".tsv"))
    #Version to deal with mouse
    elif focal_species == "Mus_musculus":
        point2 = ["GL456021", "GL456033", "GL455992", "GL455999", "GL456045", "GL456028", "JH792831", "GL456049", "GL456002", "GL456012", "GL456024", "GL456003", "KQ030493", "GL456054", "GL456053", "GL456022", "GL456017", "GL456001", "GL456060", "GL456042", "GL456026", "KB469741", "GL456044"]
        point3 = ["KB469738"]
        if l[0] in ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY", "chrM"] or l[0].endswith("_random") or l[0].startswith("chrUn"):
            out.write(com_snp.replace("REFSEQ", l[0]).replace("START", l[1]).replace("TSV_NAME", l[0] + "_" + l[1] + ".tsv"))
        elif l[0] in point2:
            out.write(com_snp.replace("REFSEQ", l[0] + ".2").replace("START", l[1]).replace("TSV_NAME", l[0] + ".2_" + l[1] + ".tsv"))
        elif l[0] in point3:
            out.write(com_snp.replace("REFSEQ", l[0] + ".3").replace("START", l[1]).replace("TSV_NAME", l[0] + ".3_" + l[1] + ".tsv"))
        else:
            out.write(com_snp.replace("REFSEQ", l[0] + ".1").replace("START", l[1]).replace("TSV_NAME", l[0] + ".1_" + l[1] + ".tsv"))
    elif focal_species == "Nycticebus_pygmaeus":
        out.write(com_snp.replace("REFSEQ", l[0]).replace("START", l[1]).replace("TSV_NAME", l[0] + "_" + l[1] + ".tsv"))
    else:
        out.write(com_snp.replace("REFSEQ", l[0] + ".1").replace("START", l[1]).replace("TSV_NAME", l[0] + ".1_" + l[1] + ".tsv"))
    if c % 4000 == 0:
        out.write("\n")
        out.write("cat *.tsv > Run" + str(tot) + ".tsv\n")
        out.write("python ../../remHead_callConvDiv.py Run" + str(tot) + ".tsv " + ignore_first + "\n")
        out.write("rm *.tsv\n")
        out.write("mv Run" + str(tot) + ".bed " + "../All")
        out.close()
        tot += 1
        os.mkdir(out_dir + "/run" + str(tot))
        out = open(out_dir + "/run" + str(tot) + "/run" + str(tot) + ".sh", 'w')
        out.write("#!/bin/bash\n#SBATCH --time=36:00:00\n#SBATCH -p hns,hbfraser\n#SBATCH --mem=8GB\n\n")

out.write("cat *.tsv > Run" + str(tot) + ".tsv\n")
out.write("python ../../remHead_callConvDiv.py Run" + str(tot) + ".tsv " + ignore_first + "\n")
out.write("rm *.tsv\n")
out.write("mv Run" + str(tot) + ".bed " + "../All\n")
out.close()
