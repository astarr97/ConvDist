import os
import sys

#Prefix of the tsv files we will pull species from
prefix = sys.argv[1]

#Lift from Species A to Species B, comma-delimited list
#Set to 0 if no liftover needed
lift = sys.argv[2]

if lift != "0":
    try:
        assert(lift.split(",")[0] in prefix)
    except:
        print("Must have the reference species in the prefix if lifting over")
        assert(False)

#Focal, Related, Outgroup and the additional species
focal = sys.argv[3]
rel = sys.argv[4]
outgroup = sys.argv[5]
add_focal = sys.argv[6]
add_rel = sys.argv[7]
add_out = sys.argv[8]

def write_beg(out):
    out.write("#!/bin/bash\n")
    out.write("#SBATCH --time=36:00:00\n")
    out.write("#SBATCH -p hbfraser\n")
    out.write("#SBATCH --mem=24GB\n")
    out.write("#SBATCH -c 1\n\n")

#Command process
command_process="python filter_and_tag_variants.py REPLACE_FILE FOCAL REL OUT ADD_FOCAL ADD_REL ADD_OUT"
command_lift="./halLiftover --bedType 3 HALFILE SRCGENOME SRCBED TGTGENOME TGTBED"
out = open("run1.sh", 'w')
write_beg(out)
c = 1
for file in os.listdir():
    if prefix in file and ".tsv" in file:
        out.write(command_process.replace("ADD_FOCAL", add_focal).replace("ADD_REL", add_rel).replace("ADD_OUT", add_out).replace("REPLACE_FILE", file).replace("FOCAL", focal).replace("REL", rel).replace("OUT", outgroup) + "\n")
        if lift != "0":
            out.write(command_lift.replace("HALFILE", "/scratch/users/astarr97/PhyloP/hg38.447way.hal").replace("SRCGENOME", lift.split(",")[0]).replace("SRCBED", file.replace(".tsv", "_ForPhyloP.bed")).replace("TGTGENOME", lift.split(",")[1]).replace("TGTBED", file.replace(".tsv", "_ForPhyloP.bed").replace(lift.split(",")[0], lift.split(",")[1])) + "\n")
        c += 1
        if c % 10 == 0:
            out.close()
            out = open("run" + str(c//10 + 1) + ".sh", 'w')
            write_beg(out)