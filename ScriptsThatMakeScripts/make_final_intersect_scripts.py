import sys
import os

#Focal species must be the first in the list!
#Usage example:
#species_list="Orcinus_orca,Bos_taurus,Sus_scrofa;Enhydra_lutris,Mustela_putorius,Mellivora_capensis;Odobenus_rosmarus,Ursus_maritimus,Spilogale_gracilis;Trichechus_manatus,Loxodonta_africana,Orycteropus_afer"
#backup_location="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/Output_Files"
#python make_final_intersect_scripts.py $species_list Aquatic do_final_intersect_Aquatic.sh $backup_location

#You can then swap everything for the Rel species with the following
#NOTE YOU MUST HAVE "Rel" in the out_file_name (argument 3)!
#species_list="Bos_taurus,Orcinus_orca,Sus_scrofa;Mustela_putorius,Enhydra_lutris,Mellivora_capensis;Ursus_maritimus,Odobenus_rosmarus,Spilogale_gracilis;Loxodonta_africana,Trichechus_manatus,Orycteropus_afer"
#backup_location="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/Output_Files"
#python make_final_intersect_scripts.py $species_list Aquatic do_final_intersect_Aquatic_Rel.sh $backup_location

species = sys.argv[1].split(";")
suffix = sys.argv[2]
out_file_name = sys.argv[3]
backup_location = sys.argv[4]

gen_folder_name = "SPECIES1_SPECIES2_NewTSVs"

intersect_initial = "bedtools intersect -sorted -wao -a SPECIES1_FILE -b SPECIES2_FILE > OUTPUT_INTERSECT"
process_initial = "python process_intersect.py OUTPUT_INTERSECT"
make_unique = "sort -u OUTPUT_READY | sort -k1,1 -k2,2n > OUTPUT_UNIQUE"

last_intersect = "bedtools intersect -sorted -wao -a OUTPUT_UNIQUE -b OUTPUT_STRANDED > OUTPUT_TO_FIX"
process_last = "python process_intersect2.py OUTPUT_TO_FIX SPECIES_LIST"

closest1 = "bedtools closest -sorted -wao -k 2 -d -a OUTPUT_TO_FIX -b SPECIES1_FILE > OUTPUT_CLOSEST1"
process_closest1 = "python process_closest1.py OUTPUT_CLOSEST1"
closest2 = "bedtools closest -sorted -wao -k 2 -d -a OUTPUT_RECLOSE -b SPECIES2_FILE > OUTPUT_CLOSEST2"
process_closest2 = "python process_closest2.py OUTPUT_CLOSEST2"

already_zipped = []
out = open(out_file_name, 'w')
c = 0
for i in range(len(species)):
    for j in range(len(species)):
        if i != j and i < j:
            c += 1
            out.write("#Doing for species pair " + str(c) + ", run block 1 for each species pair first, it will take a few hours to finish then the rest can be run quickly\n")
            out.write("#Begin block 1 for species pair " + str(c) + "\n")
            out.write('folder_name="' + gen_folder_name.replace("SPECIES1", species[i].split(",")[0]).replace("SPECIES2", species[j].split(",")[0]) + '"\n')
                       
            if "Rel" in out_file_name:
                file2 = species[j].split(",")[1] + "_" + suffix + ".FiltConv." + species[0].split(",")[1] + ".PhyloP.NearestGene.SpecSup.bed"
                if i == 0:
                    file1 = species[i].split(",")[1] + "_" + suffix + ".FiltConv.PhyloP.NearestGene.SpecSup.bed"
                else:
                    file1 = species[i].split(",")[1] + "_" + suffix + ".FiltConv." + species[0].split(",")[1] + ".PhyloP.NearestGene.SpecSup.bed"
 
                file1 = file1.replace("FiltConv.", "FiltConv.Rel.")
                print(file1)
                file2 = file2.replace("FiltConv.", "FiltConv.Rel.")
            else:
                file2 = species[j].split(",")[0] + "_" + suffix + ".FiltConv." + species[0].split(",")[0] + ".PhyloP.NearestGene.SpecSup.bed"
                if i == 0:
                    file1 = species[i].split(",")[0] + "_" + suffix + ".FiltConv.PhyloP.NearestGene.SpecSup.bed"
                else:
                    file1 = species[i].split(",")[0] + "_" + suffix + ".FiltConv." + species[0].split(",")[0] + ".PhyloP.NearestGene.SpecSup.bed"

            if "Rel" in out_file_name:
                output_intersect = species[i].split(",")[0] + "." + species[j].split(",")[0] + "." + suffix + ".Rel.Intersect.bed"
            else:
                output_intersect = species[i].split(",")[0] + "." + species[j].split(",")[0] + "." + suffix + ".Intersect.bed"
            
            out.write(intersect_initial.replace("SPECIES1_FILE", file1).replace("SPECIES2_FILE", file2).replace("OUTPUT_INTERSECT", output_intersect) + "\n")
            out.write(process_initial.replace("OUTPUT_INTERSECT", output_intersect) + "\n")
            out.write(make_unique.replace("OUTPUT_READY", output_intersect.replace("Intersect", "Ready")).replace("OUTPUT_UNIQUE", output_intersect.replace("Intersect", "Uniq")) + "\n\n")
            out.write("mkdir $folder_name\n")
            
            #We want to ignore the first column with bases if it the reference species is not in the comparison (ie if i != 0)
            if i == 0 and "Rel" not in out_file_name:
                ignore_first = "0"
                species_list = species[j].split(",")[0] + "," + ",".join(species[i].split(",")[1:]) + "," + ",".join(species[j].split(",")[1:])
            elif "Rel" not in out_file_name:
                ignore_first = "1"
                species_list = species[i].split(",")[0] + "," + species[j].split(",")[0] + "," + ",".join(species[i].split(",")[1:]) + "," + ",".join(species[j].split(",")[1:])
            elif "Rel" in out_file_name and i == 0:
                ignore_first = "2"
                species_list = species[i].split(",")[0] + "," + species[j].split(",")[0] + "," + ",".join(species[i].split(",")[2:]) + "," + ",".join(species[j].split(",")[1:])
            elif "Rel" in out_file_name:
                ignore_first = "1"
                species_list = species[i].split(",")[0] + "," + species[j].split(",")[0] + "," + ",".join(species[i].split(",")[1:]) + "," + ",".join(species[j].split(",")[1:])
            else:
                print("Must be an error")
            if "Rel" in out_file_name:
                out.write("python make_new_tsv_scripts.py " + output_intersect.replace("Intersect", "Uniq") + " " + species[0].split(",")[1] + " " + species_list + " $folder_name " + ignore_first + "\n")
            else:
                out.write("python make_new_tsv_scripts.py " + output_intersect.replace("Intersect", "Uniq") + " " + species[0].split(",")[0] + " " + species_list + " $folder_name " + ignore_first + "\n")
            out.write("cd $folder_name\ncp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh ./\n./driver.sh\ncd ..\n\n")
            
            out.write("#Begin block 2 for species pair " + str(c) + "\n")
            out.write("cat ${folder_name}/All/*.bed | sort -k1,1 -k2,2n > " + output_intersect.replace("Intersect", "Stranded") + "\n")
            
            out.write(last_intersect.replace("OUTPUT_UNIQUE", output_intersect.replace("Intersect", "Uniq")).replace("OUTPUT_STRANDED", output_intersect.replace("Intersect", "Stranded")).replace("OUTPUT_TO_FIX", output_intersect.replace("Intersect", "ToFix")) + "\n")
            
            species_list_pass = species[i].split(",")[0] + "," + species[j].split(",")[0] + "," + ",".join(species[i].split(",")[1:]) + "," + ",".join(species[j].split(",")[1:])
            out.write(process_last.replace("OUTPUT_TO_FIX", output_intersect.replace("Intersect", "ToFix")).replace("SPECIES_LIST", species_list_pass) + "\n\n")
            
            #Now we want to do the closest site approach
            out.write("#Begin block 3 for species pair " + str(c) + "\n")
            out.write(closest1.replace("OUTPUT_TO_FIX", output_intersect.replace("Intersect", "ToFix")).replace("SPECIES1_FILE", file1).replace("OUTPUT_CLOSEST1", output_intersect.replace("Intersect", species[i].split(",")[0] + "Closest")) + "\n")
            out.write(process_closest1.replace("OUTPUT_CLOSEST1", output_intersect.replace("Intersect", species[i].split(",")[0] + "Closest")) + "\n")
            out.write(closest2.replace("OUTPUT_RECLOSE", output_intersect.replace("Intersect", species[i].split(",")[0] + "ReClose")).replace("SPECIES2_FILE", file2).replace("OUTPUT_CLOSEST2", output_intersect.replace("Intersect", species[i].split(",")[0] + "Closest." + species[j].split(",")[0] + "Closest")) + "\n")
            out.write(process_closest2.replace("OUTPUT_CLOSEST2", output_intersect.replace("Intersect", species[i].split(",")[0] + "Closest." + species[j].split(",")[0] + "Closest")) + "\n")
            
            if "Rel" in out_file_name:
                final_file_closest = species[i].split(",")[1] + "." + species[j].split(",")[1] + "." + suffix + ".Rel.ClosestVar.bed"
            else:
                final_file_closest = species[i].split(",")[0] + "." + species[j].split(",")[0] + "." + suffix + ".ClosestVar.bed"
            out.write("mv " + output_intersect.replace("Intersect", species[i].split(",")[0] + "." + species[j].split(",")[0] + "ClosReady") + " " + final_file_closest + "\n\n")
            
            out.write("#Begin block 4 for species pair " + str(c) + "\n")
            #Now we want to compress and back everything up
            out.write("gzip -c " + final_file_closest + " > " + final_file_closest + ".gz\n")
            if "Rel" in out_file_name:
                out.write("mv " + output_intersect.replace(".Intersect.bed", ".ConvDiv.txt") + " " + output_intersect.replace(".Intersect.bed", ".ConvDiv.txt").replace(species[i].split(",")[0], species[i].split(",")[1]).replace(species[j].split(",")[0], species[j].split(",")[1]) + "\n")
                out.write("gzip -c " + output_intersect.replace(".Intersect.bed", ".ConvDiv.txt").replace(species[i].split(",")[0], species[i].split(",")[1]).replace(species[j].split(",")[0], species[j].split(",")[1]) + " > " + output_intersect.replace(".Intersect.bed", ".ConvDiv.txt").replace(species[i].split(",")[0], species[i].split(",")[1]).replace(species[j].split(",")[0], species[j].split(",")[1]) + ".gz\n")
            else:
                out.write("gzip -c " + output_intersect.replace(".Intersect.bed", ".ConvDiv.txt") + " > " + output_intersect.replace(".Intersect.bed", ".ConvDiv.txt") + ".gz\n")
            if file1 not in already_zipped:
                out.write("gzip -c " + file1 + " > " + file1 + ".gz\n")
                out.write("gzip -c " + file1.replace("FiltConv", "FiltPoly") + " > " + file1.replace("FiltConv", "FiltPoly") + ".gz\n")
                already_zipped.append(file1)
            if file2 not in already_zipped:
                out.write("gzip -c " + file2 + " > " + file2 + ".gz\n")
                out.write("gzip -c " + file2.replace("FiltConv", "FiltPoly") + " > " + file2.replace("FiltConv", "FiltPoly") + ".gz\n")
                already_zipped.append(file2)
            out.write("\n\n")
try:
    os.listdir(backup_location)
except:
    print("Backup location does not exist")
    assert(False)
out.write("mv *.gz " + backup_location)
out.close()

            