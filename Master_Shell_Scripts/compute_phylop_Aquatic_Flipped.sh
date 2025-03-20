ml load biology samtools
ml load python

to_mask="Tragulus_javanicus,Moschus_moschiferus,Bubalus_bubalis,Bos_taurus,Bos_indicus,Bos_mutus,Bison_bison,Beatragus_hunteri,Ammotragus_lervia,Hemitragus_hylocrius,Capra_hircus,Capra_aegagrus,Ovis_aries,Ovis_canadensis,Pantholops_hodgsonii,Saiga_tatarica,Okapia_johnstoni,Giraffa_tippelskirchi,Antilocapra_americana,Odocoileus_virginianus,Rangifer_tarandus,Elaphurus_davidianus,Spilogale_gracilis,Mustela_putorius,Ailurus_fulgens,Mellivora_capensis,Enhydra_lutris,Pteronura_brasiliensis,Loxodonta_africana"
phylop_species="Orcinus_orca"
folder_name="Orcinus_orca_PhyloP_MaskAquaticFlipped"
#cp -r /home/groups/hbfraser/Common_Software /scratch/users/astarr97/Common_Software
#mkdir /scratch/users/astarr97/astarr_scripts
#cp -r /home/groups/hbfraser/astarr_scripts/AccelConvDist /scratch/users/astarr97/astarr_scripts/AccelConvDist

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta /scratch/users/astarr97/PhyloP/hg38.447way.hal $phylop_species > $phylop_species.fasta
samtools faidx $phylop_species.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $phylop_species.fasta.fai | sort -n -r -k 3,3 > ${phylop_species}_contigs.bed

python make_scripts_mask_phylop.py ${phylop_species}_contigs.bed $to_mask $phylop_species $folder_name

cp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh /scratch/users/astarr97/PhyloP/Feasibility/$folder_name
cd $folder_name
./driver.sh
cd ..

#Catted all the resulting bed files together, sorted, then summed by chromosome to validate that we actually cover the entire orca genome.
#python check_all.py
#sort -k1,1 -k2,2n Orcinus_orca_PhyloP_MaskAquatic_New_All.bed > Orcinus_orca_PhyloP_MaskAquatic_New_All.sort.bed
