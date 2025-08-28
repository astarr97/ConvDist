import sys
import copy

#The first argument is the config file, the second argument is the reference species, and the third argument is whether this is being done for the focal species (argument should be 0) or rel (1)
config = sys.argv[1]
ref_species = sys.argv[2]
rel = int(sys.argv[3])

trait = config.replace("Config_GetDL_", "").replace(".csv", "")
print(trait)
o = open(config)
if rel:
    out1 = open("script_Rel_" + trait + "_1.sh", 'w')
    out2 = open("script_Rel_" + trait + "_2.sh", 'w')
    out3 = open("script_Rel_" + trait + "_3.sh", 'w')
else:
    out1 = open("script_" + trait + "_1.sh", 'w')
    out2 = open("script_" + trait + "_2.sh", 'w')
    out3 = open("script_" + trait + "_3.sh", 'w')
    
out0 = open("script_" + trait + "_0.sh", 'w')
#Write out headers
out1.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=16GB\n\n")
out3.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=16GB\n\n")

#This is a string that encodes all the relevant variables at the beginning of the scripts
var_string = "#This is the reference species\nref_species=REFSPEC\n#This is the first focal species\nfoc1=FOCAL1\n#And the second focal species\nfoc2=FOCAL2\n#This is the first rel species\nrel1=REL1\n#And the second rel species\nrel2=REL2\n#This is the first outgroup species\nout1=OUT1\n#This is the second outgroup species\nout2=OUT2\n#This is the LCA for the first focal species clade\nfoc1_lca=FOC1_LCA\n#This is the LCA for the second focal species clade\nfoc2_lca=FOC2_LCA\n#This is the trait of interest\ntrait=TRAIT\n#This is the folder on OAK to back up to\nbackup=/oak/stanford/groups/hbfraser/astarr/AccelConvDist/${trait}/Output_Files_DL\n\n"

#Command to make fasta file and contig files
com_get_fasta = "/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta /scratch/users/astarr97/PhyloP/hg38.447way.hal GENOME > GENOME.fasta\nsamtools faidx GENOME.fasta\nawk 'BEGIN {FS=" + '"\\t"' + "}; {print $1 FS " + '"0"' + " FS $2}' GENOME.fasta.fai | sort -n -r -k 3,3 > GENOME_contigs.bed\n\n"

#Command to make bed files from the ClosestVar.bed.gz files
com_make_beds = "python make_beds_2114.py FOC1.FOC2.TRAIT.ClosestVar.bed.gz REF_SPECIES RELLL\n\n"

#Command to do halLiftover
com_lift = "/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal REF_SPECIES FILE NEW_GENOME OUT_ILE"

#Command to make new TSVs for each of the four bed files
com_tsv1="folder=${foc1_lca}_${foc1}_${foc2}_NewTSVs\nmkdir $folder\npython make_new_tsv_scripts_dl.py $foc1.$foc2.$trait.Focal.${foc1_lca}.bed ${foc1_lca} ${foc1},${foc2},${rel1},${out1},${rel2},${out2} $folder\ncd $folder\ncp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh ./\n./driver.sh\ncd ..\n\n"
com_tsv2="folder=${foc2_lca}_${foc1}_${foc2}_NewTSVs\nmkdir $folder\npython make_new_tsv_scripts_dl.py $foc1.$foc2.$trait.Focal.${foc2_lca}.bed ${foc2_lca} ${foc1},${foc2},${rel1},${out1},${rel2},${out2} $folder\ncd $folder\ncp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh ./\n./driver.sh\ncd ..\n\n"
com_tsv3="folder=${foc1_lca}_${foc1}_${foc2}_Closest_NewTSVs\nmkdir $folder\npython make_new_tsv_scripts_dl.py $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.bed ${foc1_lca} ${foc1},${foc2},${rel1},${out1},${rel2},${out2} $folder\ncd $folder\ncp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh ./\n./driver.sh\ncd ..\n\n"
com_tsv4="folder=${foc2_lca}_${foc1}_${foc2}_Closest_NewTSVs\nmkdir $folder\npython make_new_tsv_scripts_dl.py $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.bed ${foc2_lca} ${foc1},${foc2},${rel1},${out1},${rel2},${out2} $folder\ncd $folder\ncp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh ./\n./driver.sh\ncd ..\n\n"

seen = []

for line in o:
    l = line.replace("\n", "").split(",")
    var_string2 = copy.deepcopy(var_string)
    var_string2 = var_string2.replace("REFSPEC", ref_species).replace("FOCAL1", l[1]).replace("FOCAL2", l[5]).replace("REL1", l[2]).replace("REL2", l[6]).replace("OUT1", l[3]).replace("OUT2", l[7]).replace("FOC1_LCA", l[0]).replace("FOC2_LCA", l[4]).replace("TRAIT", trait)
    out1.write(var_string2)
    out2.write(var_string2)
    out3.write(var_string2)
    
    
    if l[0] not in seen:
        out0.write("trait=" + trait + "\n")
        out0.write("backup=/oak/stanford/groups/hbfraser/astarr/AccelConvDist/${trait}/Output_Files_DL\n\n")
        out0.write("mkdir $backup\n\n")
        out0.write("foc1_lca=" + l[0] + "\n")
        seen.append(l[0])
        out0.write(com_get_fasta.replace("GENOME", "${foc1_lca}"))
    if l[4] not in seen:
        out0.write("foc2_lca=" + l[4] + "\n")
        seen.append(l[4])
        out0.write(com_get_fasta.replace("GENOME", "${foc2_lca}"))
    
    
    
    if rel:
        out1.write(com_make_beds.replace("FOC1", "$foc1").replace("FOC2", "$foc2").replace("TRAIT", "$trait").replace("REF_SPECIES", ref_species).replace("RELLL", "1").replace("ClosestVar", "Rel.ClosestVar"))
        out1.write("#Lift over each file to the respective ancestral genomes\n")
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Rel.Focal.bed").replace("NEW_GENOME", l[0]).replace("OUT_ILE", "$foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.bed\n"))
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Rel.Focal.bed").replace("NEW_GENOME", l[4]).replace("OUT_ILE", "$foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.bed\n"))
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Rel.Closest_S1.bed").replace("NEW_GENOME", l[0]).replace("OUT_ILE", "$foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.bed\n"))
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Rel.Closest_S2.bed").replace("NEW_GENOME", l[4]).replace("OUT_ILE", "$foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.bed\n") + "\n\n")

    else:
        out1.write(com_make_beds.replace("FOC1", "$foc1").replace("FOC2", "$foc2").replace("TRAIT", "$trait").replace("REF_SPECIES", ref_species).replace("RELLL", "0"))
        out1.write("#Lift over each file to the respective ancestral genomes\n")
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Focal.bed").replace("NEW_GENOME", l[0]).replace("OUT_ILE", "$foc1.$foc2.$trait.Focal.${foc1_lca}.bed\n"))
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Focal.bed").replace("NEW_GENOME", l[4]).replace("OUT_ILE", "$foc1.$foc2.$trait.Focal.${foc2_lca}.bed\n"))
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Closest_S1.bed").replace("NEW_GENOME", l[0]).replace("OUT_ILE", "$foc1.$foc2.$trait.Closest_S1.${foc1_lca}.bed\n"))
        out1.write(com_lift.replace("REF_SPECIES", "$ref_species").replace("FILE", "$foc1.$foc2.$trait.Closest_S2.bed").replace("NEW_GENOME", l[4]).replace("OUT_ILE", "$foc1.$foc2.$trait.Closest_S2.${foc2_lca}.bed\n") + "\n\n")
    
    if rel:
        out2.write(com_tsv1.replace("${foc1},${foc2},${rel1},${out1},${rel2},${out2}", "${rel1},${rel2},${foc1},${out1},${foc2},${out2}").replace("$foc1.$foc2.$trait.", "$foc1.$foc2.$trait.Rel.").replace("_NewTSVs", "_Rel_NewTSVs"))
        out2.write(com_tsv2.replace("${foc1},${foc2},${rel1},${out1},${rel2},${out2}", "${rel1},${rel2},${foc1},${out1},${foc2},${out2}").replace("$foc1.$foc2.$trait.", "$foc1.$foc2.$trait.Rel.").replace("_NewTSVs", "_Rel_NewTSVs"))
        out2.write(com_tsv3.replace("${foc1},${foc2},${rel1},${out1},${rel2},${out2}", "${rel1},${rel2},${foc1},${out1},${foc2},${out2}").replace("$foc1.$foc2.$trait.", "$foc1.$foc2.$trait.Rel.").replace("_NewTSVs", "_Rel_NewTSVs"))
        out2.write(com_tsv4.replace("${foc1},${foc2},${rel1},${out1},${rel2},${out2}", "${rel1},${rel2},${foc1},${out1},${foc2},${out2}").replace("$foc1.$foc2.$trait.", "$foc1.$foc2.$trait.Rel.").replace("_NewTSVs", "_Rel_NewTSVs"))
    else:
        out2.write(com_tsv1)
        out2.write(com_tsv2)
        out2.write(com_tsv3)
        out2.write(com_tsv4)
    
    if rel:
        #Commands for the third script
        out3.write("folder=${foc1_lca}_${foc1}_${foc2}_Rel_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.bed All_${folder}.tsv ${rel1},${rel2},${foc1},${foc2},${foc1_lca}\n\n")
        
        out3.write("folder=${foc2_lca}_${foc1}_${foc2}_Rel_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.bed All_${folder}.tsv ${rel1},${rel2},${foc1},${foc2},${foc2_lca}\n\n")
        
        out3.write("folder=${foc1_lca}_${foc1}_${foc2}_Closest_Rel_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.bed All_${folder}.tsv ${rel1},${foc1},${foc1_lca}\n\n")
        
        out3.write("folder=${foc2_lca}_${foc1}_${foc2}_Closest_Rel_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.bed All_${folder}.tsv ${rel2},${foc2},${foc2_lca}\n\n")
        
        out3.write("#Expand to 2114 bases on either side and get fasta for focal\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.Comped.bed ${foc1_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc1_lca}.fasta -bed $foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.Comped.2114.fastabed\n\n")
        
        out3.write("#Do the same for other species for focal\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.Comped.bed ${foc2_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc2_lca}.fasta -bed $foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.Comped.2114.fastabed\n\n")
        
        out3.write("#Do the same for Closest_S1\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.Comped.bed ${foc1_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc1_lca}.fasta -bed $foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.Comped.2114.fastabed\n\n")
        
        out3.write("#Do the same for Closest_S2\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.Comped.bed ${foc2_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc2_lca}.fasta -bed $foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.Comped.2114.fastabed\n\n")
        
        #Now generate all the input files, compress, and backup
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.Comped.2114.fastabed $foc1\n")
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.Comped.2114.fastabed $foc2\n")
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.Comped.2114.fastabed $foc1\n")
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.Comped.2114.fastabed $foc2\n\n")
        
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.$foc1.Final.fastabed > $backup/$foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.$foc1.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.$foc1.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Rel.Focal.${foc1_lca}.$foc1.ToCompare.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.$foc2.Final.fastabed > $backup/$foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.$foc2.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.$foc2.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Rel.Focal.${foc2_lca}.$foc2.ToCompare.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.$foc1.Final.fastabed > $backup/$foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.$foc1.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.$foc1.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Rel.Closest_S1.${foc1_lca}.$foc1.ToCompare.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.$foc2.Final.fastabed > $backup/$foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.$foc2.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.$foc2.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Rel.Closest_S2.${foc2_lca}.$foc2.ToCompare.fastabed.gz\n\n")
    else:
        out3.write("folder=${foc1_lca}_${foc1}_${foc2}_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Focal.${foc1_lca}.bed All_${folder}.tsv ${foc1},${foc2},${rel1},${rel2},${foc1_lca}\n\n")
        
        out3.write("folder=${foc2_lca}_${foc1}_${foc2}_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Focal.${foc2_lca}.bed All_${folder}.tsv ${foc1},${foc2},${rel1},${rel2},${foc2_lca}\n\n")
        
        out3.write("folder=${foc1_lca}_${foc1}_${foc2}_Closest_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.bed All_${folder}.tsv ${foc1},${rel1},${foc1_lca}\n\n")
        
        out3.write("folder=${foc2_lca}_${foc1}_${foc2}_Closest_NewTSVs\n")
        out3.write("cat $folder/All/*tsv > All_${folder}.tsv\n")
        out3.write("python correct_strand.py $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.bed All_${folder}.tsv ${foc2},${rel2},${foc2_lca}\n\n")
        
        out3.write("#Expand to 2114 bases on either side and get fasta for focal\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Focal.${foc1_lca}.Comped.bed ${foc1_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc1_lca}.fasta -bed $foc1.$foc2.$trait.Focal.${foc1_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Focal.${foc1_lca}.Comped.2114.fastabed\n\n")
        
        out3.write("#Do the same for other species for focal\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Focal.${foc2_lca}.Comped.bed ${foc2_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc2_lca}.fasta -bed $foc1.$foc2.$trait.Focal.${foc2_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Focal.${foc2_lca}.Comped.2114.fastabed\n\n")
        
        out3.write("#Do the same for Closest_S1\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.Comped.bed ${foc1_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc1_lca}.fasta -bed $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.Comped.2114.fastabed\n\n")
        
        out3.write("#Do the same for Closest_S2\n")
        out3.write("python expand_2114.py $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.Comped.bed ${foc2_lca}_contigs.bed\n")
        out3.write("bedtools getfasta -bedOut -name -fi ${foc2_lca}.fasta -bed $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.Comped.2114.bed > $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.Comped.2114.fastabed\n\n")
        
        #Now generate all the input files, compress, and backup
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Focal.${foc1_lca}.Comped.2114.fastabed $foc1\n")
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Focal.${foc2_lca}.Comped.2114.fastabed $foc2\n")
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.Comped.2114.fastabed $foc1\n")
        out3.write("python gen_final_input.py $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.Comped.2114.fastabed $foc2\n\n")
        
        out3.write("gzip -c $foc1.$foc2.$trait.Focal.${foc1_lca}.$foc1.Final.fastabed > $backup/$foc1.$foc2.$trait.Focal.${foc1_lca}.$foc1.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Focal.${foc1_lca}.$foc1.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Focal.${foc1_lca}.$foc1.ToCompare.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Focal.${foc2_lca}.$foc2.Final.fastabed > $backup/$foc1.$foc2.$trait.Focal.${foc2_lca}.$foc2.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Focal.${foc2_lca}.$foc2.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Focal.${foc2_lca}.$foc2.ToCompare.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.$foc1.Final.fastabed > $backup/$foc1.$foc2.$trait.Closest_S1.${foc1_lca}.$foc1.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Closest_S1.${foc1_lca}.$foc1.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Closest_S1.${foc1_lca}.$foc1.ToCompare.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.$foc2.Final.fastabed > $backup/$foc1.$foc2.$trait.Closest_S2.${foc2_lca}.$foc2.Final.fastabed.gz\n")
        out3.write("gzip -c $foc1.$foc2.$trait.Closest_S2.${foc2_lca}.$foc2.ToCompare.fastabed > $backup/$foc1.$foc2.$trait.Closest_S2.${foc2_lca}.$foc2.ToCompare.fastabed.gz\n\n")

out1.close()
out2.close()
out3.close()

