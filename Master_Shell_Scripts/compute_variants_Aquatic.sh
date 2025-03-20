focal_species="Orcinus_orca"

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta /scratch/users/astarr97/PhyloP/hg38.447way.hal $focal_species > $focal_species.fasta
samtools faidx $focal_species.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $focal_species.fasta.fai | sort -n -r -k 3,3 > ${focal_species}_contigs.bed

focal_species="Odobenus_rosmarus"

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta /scratch/users/astarr97/PhyloP/hg38.447way.hal $focal_species > $focal_species.fasta
samtools faidx $focal_species.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $focal_species.fasta.fai | sort -n -r -k 3,3 > ${focal_species}_contigs.bed

focal_species="Enhydra_lutris"

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta /scratch/users/astarr97/PhyloP/hg38.447way.hal $focal_species > $focal_species.fasta
samtools faidx $focal_species.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $focal_species.fasta.fai | sort -n -r -k 3,3 > ${focal_species}_contigs.bed

focal_species="Trichechus_manatus"

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta /scratch/users/astarr97/PhyloP/hg38.447way.hal $focal_species > $focal_species.fasta
samtools faidx $focal_species.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $focal_species.fasta.fai | sort -n -r -k 3,3 > ${focal_species}_contigs.bed

python make_scripts_get_variants_dist.py Config_Aquatic.txt Aquatic

folder_name="Orcinus_orca_Aquatic"
cp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh /scratch/users/astarr97/PhyloP/Feasibility/$folder_name
cd $folder_name
./driver.sh
cd ..

folder_name="Odobenus_rosmarus_Aquatic"
cp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh /scratch/users/astarr97/PhyloP/Feasibility/$folder_name
cd $folder_name
./driver.sh
cd ..

folder_name="Enhydra_lutris_Aquatic"
cp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh /scratch/users/astarr97/PhyloP/Feasibility/$folder_name
cd $folder_name
./driver.sh
cd ..

folder_name="Trichechus_manatus_Aquatic"
cp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh /scratch/users/astarr97/PhyloP/Feasibility/$folder_name
cd $folder_name
./driver.sh
cd ..