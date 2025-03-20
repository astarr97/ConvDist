import sys
import os
import pandas as pd

#This is a configuration file, each row should generally correspond to an independent evolution of the trait
#The first column is the set of species with the trait
#The second columns is the set of closely related species without the trait
#The third columns is the set of outgroup species
#The fourth column is the bed file of interest

#In each of the first three columns, all species after the semi-colon will be pulled but will not be included when computing things during later steps!
config_file = sys.argv[1]

#This is the suffix of the folder to output everything to
folder_suffix = sys.argv[2]

#Prefix for the resulting files
out_prefix = "SNPS_"

config = pd.read_csv(config_file, sep = "\t")
print(config.columns)

#Function to write script headers
def write_beg(out):
    out.write("#!/bin/bash\n")
    out.write("#SBATCH --time=72:00:00\n")
    out.write("#SBATCH -p hbfraser\n")
    out.write("#SBATCH --mem=24GB\n")
    out.write("#SBATCH -c 1\n\n")
    
def write_block(out, refseq, s_end, start, length, fout):
    rel_bed = fout.replace(".tsv", "_" + spec_rel + ".bed")
    foc_bed = fout.replace(".tsv", "_" + spec_focus + ".bed")
    if length != 2000000:
        out.write(com_snp.replace("REFSEQ", refseq).replace("START", str(start)).replace("LENGTH", str(length - 1)).replace("OUT_FILE", fout))
    else:
        out.write(com_snp.replace("REFSEQ", refseq).replace("START", str(start)).replace("LENGTH", str(length)).replace("OUT_FILE", fout))
    out.write("\n")

#Function to write out the commands to run at the end of a script
#spec_rel is unused
def write_end(out, run_num, spec_rel):
    out_tsv = "Run" + str(run_num) + "_ALL.tsv"
    out.write("cat *.tsv > " + out_tsv + "\n")
    #First one is for the focal species
    out.write(com_filt_snp.replace("FILE", out_tsv).replace("SPEC_FOCUS", spec_focus).replace("SPEC_REL", spec_rel).replace("SPEC_OUT", spec_out))
    
    #Ignore
    #out.write(com_filt_snp2.replace("FILE", out_tsv).replace("SPEC_FOCUS", spec_focus).replace("SPEC_REL", spec_rel).replace("SPEC_OUT", spec_out))
    out.write("mv " + out_tsv + " ../All\n")
    out.write("mv " + out_tsv.replace(".tsv", "_ForPhyloP.bed") + " ../All\n")
    #Ignore
    #out.write("mv " + out_tsv.replace(".tsv", "_" + spec_rel + "_ForPhyloP.bed" + " /All_" + spec_rel))
    out.write("\n")

def check_commas(s):
    if not s:
        return s
    s = s.replace(",,,", ",")
    s = s.replace(",,", ",")
    if s[0] == ",":
        s = s[1:]
        if not s:
            return s
    if s[-1] == ",":
        s = s[:-1]
        if not s:
            return s
    return s
    
for index, row in config.iterrows():
    #spec_focus is the species that has the trait of interest
    spec_focus = row["Focal_species"].split(",")[0]
    if len(spec_focus.split(";")) != 1:
        spec_focus = spec_focus.split(";")[0]
    #spec_rel is the closest living relative of the spec_focus that does not have the trait of interest
    spec_rel = row["Related_species"].split(",")[0]
    if len(spec_rel.split(";")) != 1:
        spec_rel = spec_rel.split(";")[0]
    #out_species is the out_group used to determine if a trait is spec_focus-derived or comp_species-derived
    spec_out = row["Outgroup_species"].split(",")[0]
    if len(spec_out.split(";")) != 1:
        spec_out = spec_out.split(";")[0]
    print(spec_focus, spec_rel, spec_out)
    
    if ";" in row["Focal_species"]:
        spec_add_focus = ",".join(row["Focal_species"].split(";")[0].split(",")[1:]) + "," + row["Focal_species"].split(";")[1]
        spec_focus_ignore = row["Focal_species"].split(";")[1]
        spec_focus_use = ",".join(row["Focal_species"].split(";")[0].split(",")[1:])
    else:
        spec_add_focus = ",".join(row["Focal_species"].split(",")[1:])
        spec_focus_use = ",".join(row["Focal_species"].split(",")[1:])
        spec_focus_ignore = "0"
    if ";" in row["Related_species"]:
        spec_add_rel = ",".join(row["Related_species"].split(";")[0].split(",")[1:]) + "," + row["Related_species"].split(";")[1]
        spec_rel_ignore = row["Related_species"].split(";")[1]
        spec_rel_use = ",".join(row["Related_species"].split(";")[0].split(",")[1:])
    else:
        spec_add_rel = ",".join(row["Related_species"].split(",")[1:])
        spec_rel_use = ",".join(row["Related_species"].split(",")[1:])
        spec_rel_ignore = "0"
    if ";" in row["Outgroup_species"]:
        spec_add_out = ",".join(row["Outgroup_species"].split(";")[0].split(",")[1:]) + "," + row["Outgroup_species"].split(";")[1]
        spec_out_ignore = row["Outgroup_species"].split(";")[1]
        spec_out_use = ",".join(row["Outgroup_species"].split(";")[0].split(",")[1:])
    else:
        spec_add_out = ",".join(row["Outgroup_species"].split(",")[1:])
        spec_out_use = ",".join(row["Outgroup_species"].split(",")[1:])
        spec_out_ignore = "0"
    
    #Check to make sure that we have removed all extraneous commas
    spec_add_focus = check_commas(spec_add_focus)
    spec_add_rel = check_commas(spec_add_rel)
    spec_add_out = check_commas(spec_add_out)
    spec_add = spec_add_focus + "," + spec_add_rel + "," + spec_add_out
    spec_add = check_commas(spec_add)
    
    spec_focus_use = check_commas(spec_focus_use)
    spec_rel_use = check_commas(spec_rel_use)
    spec_out_use = check_commas(spec_out_use)
    spec_focus_ignore = check_commas(spec_focus_ignore)
    spec_rel_ignore = check_commas(spec_rel_ignore)
    spec_out_ignore = check_commas(spec_out_ignore)
    
    print(spec_add)
    #Command to pull SNPs
    com_snp = "/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halSnps --minSpeciesForSnp 1 --refSequence REFSEQ --start START --length LENGTH --tsv OUT_FILE /scratch/users/astarr97/PhyloP/hg38.447way.hal " + spec_focus + " " + ",".join([spec_rel, spec_out]) + "," + spec_add + "\n"
    
    #Command to filter SNPs
    ### NEED TO ADJUST TO PASS ARGUMENT FOR WHICH SPECIES TO IGNORE! ###
    #, spec_focus_ignore, spec_rel_ignore, spec_out_ignore don't think we need this for now, can copy and paste after "spec_out_use" in below line if we do end up needing it
    if not spec_focus_use:
        spec_focus_use = "0"
    if not spec_rel_use:
        spec_rel_use = "0"
    if not spec_out_use:
        spec_out_use = "0"
        
    com_filt_snp = "python /scratch/users/astarr97/astarr_scripts/AccelConvDist/filter_and_tag_variants.py FILE SPEC_FOCUS SPEC_REL SPEC_OUT " + " ".join([spec_focus_use, spec_rel_use, spec_out_use]) + "\n"
    #Ignore
    #com_filt_snp2 = "python /scratch/users/astarr97/astarr_scripts/AccelConvDist/filter_and_tag_variants.py FILE SPEC_REL SPEC_FOCUS SPEC_OUT " + " ".join([spec_rel_use, spec_focus_use, spec_out_use]) + " 1\n"
    os.mkdir(spec_focus + "_" + folder_suffix)
    os.mkdir(spec_focus + "_" + folder_suffix + "/run1")
    os.mkdir(spec_focus + "_" + folder_suffix + "/All")
    #Ignore
    #os.mkdir(spec_focus + "_" + folder_suffix + "/All_" + spec_rel)
    out = open(spec_focus + "_" + folder_suffix + "/run1/" + "run1.sh", 'w')
    write_beg(out)
    base_sum = 0
    c = 1
    
    o = open(row["Contigs_file"])
    run_contigs = 0
    for line in o:
        run_contigs += 1
        l = line.replace("\n", "").split("\t")
        cur_len = int(l[2])
        #Iterate through in batches of 2,000,000 bases, resetting and opening a new script if we reach 10,000,000 bases
        if cur_len < 2000000:
            fout = out_prefix + l[0] + "_" + str(0) + "-" + str(int(l[2])) + ".tsv"
            write_block(out, l[0], int(l[2]), 0, int(l[2]), fout)
            base_sum += cur_len
        else:
            s_end = 2000000
            while s_end < cur_len:
                if base_sum > 10000000:
                    write_end(out, c, spec_rel)
                    out.close()
                    c += 1
                    os.mkdir(spec_focus + "_" + folder_suffix + "/run" + str(c))
                    out = open(spec_focus + "_" + folder_suffix + "/run" + str(c) + "/" + "run" + str(c) + ".sh", 'w')
                    write_beg(out)
                    base_sum = 0
                    run_contigs = 0
                fout = out_prefix + l[0] + "_" + str(s_end - 2000000) + "-" + str(s_end) + ".tsv"
                write_block(out, l[0], s_end, s_end - 2000000, 2000000, fout)
                base_sum += 2000000
                s_end += 2000000
            fout = out_prefix + l[0] + "_" + str(s_end - 2000000) + "-" + str(l[2]) + ".tsv"
            write_block(out, l[0], s_end, s_end - 2000000, int(l[2]) - (s_end - 2000000), fout)
            base_sum += int(l[2]) - (s_end - 2000000)
            
    
        if base_sum > 10000000 or run_contigs >= 750:
            write_end(out, c, spec_rel)
            out.close()
            c += 1
            os.mkdir(spec_focus + "_" + folder_suffix + "/run" + str(c))
            out = open(spec_focus + "_" + folder_suffix + "/run" + str(c) + "/" + "run" + str(c) + ".sh", 'w')
            write_beg(out)
            base_sum = 0
            run_contigs = 0
    write_end(out, c, spec_rel)
    out.close()
