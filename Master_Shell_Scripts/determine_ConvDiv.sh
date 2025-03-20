cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/ScriptsThatMakeScripts/make_final_intersect_scripts.py ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/ScriptsThatMakeScripts/make_new_tsv_scripts.py ./

cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/process_closest1.py ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/process_closest2.py ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/process_intersect.py ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/process_intersect2.py ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/remHead_callConvDiv.py ./

#Focal species must be the first in the list!
#PhyloP reference species must be the very first species (second species for rel)!
#Usage example:
species_list="Orcinus_orca,Bos_taurus,Sus_scrofa;Enhydra_lutris,Mustela_putorius,Mellivora_capensis;Odobenus_rosmarus,Ursus_maritimus,Spilogale_gracilis;Trichechus_manatus,Loxodonta_africana,Orycteropus_afer"
backup_location="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/Output_Files"
python make_final_intersect_scripts.py $species_list Aquatic do_final_intersect_Aquatic.sh $backup_location

#You can then swap everything for the Rel species with the following
#NOTE YOU MUST HAVE "Rel" in the out_file_name (argument 3)!
species_list="Bos_taurus,Orcinus_orca,Sus_scrofa;Mustela_putorius,Enhydra_lutris,Mellivora_capensis;Spilogale_gracilis,Odobenus_rosmarus,Ursus_maritimus;Loxodonta_africana,Trichechus_manatus,Orycteropus_afer"
backup_location="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/Output_Files"
python make_final_intersect_scripts.py $species_list Aquatic do_final_intersect_Aquatic_Rel.sh $backup_location

### THIS DOES NOT CURRENTLY WORK, WORKING ON FIXING IT ###
### NOTE IF THE REFERENCE PHYLOP SPECIES IS NOT ONE OF THE FOCAL SPECIES, YOU MUST INCLUDE A BURN-IN AS BELOW ###
#The Mus_musculus acts as a burn in since everything is referenced to Mus_musculus, not Acomys_cahirinus
#This is a somewhat hacky way to get the script to realize that none of the focal species are the reference species and that Mus_musculus is the reference species
species_list="Mus_musculus,Mus_musculus,Mus_musculus;Acomys_cahirinus,Psammomys_obesus,Mus_musculus;Dasypus_novemcinctus,Tamandua_tetradactyla,Loxodonta_africana"
backup_location="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Armor/Output_Files"
python make_final_intersect_scripts.py $species_list Armor do_final_intersect_Armor.sh $backup_location

#You can do the same thing with Rel.
#You can then swap everything for the Rel species with the following
#NOTE YOU MUST HAVE "Rel" in the out_file_name (argument 3)!
species_list="Mus_musculus,Mus_musculus,Mus_musculus;Psammomys_obesus,Acomys_cahirinus,Mus_musculus;Tamandua_tetradactyla,Dasypus_novemcinctus,Loxodonta_africana"
backup_location="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Armor/Output_Files"
python make_final_intersect_scripts.py $species_list Armor do_final_intersect_Armor_Rel.sh $backup_location
