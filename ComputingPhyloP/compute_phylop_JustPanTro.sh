to_mask="Pan_troglodytes,Pan_paniscus"

phylop_species="Pan_troglodytes"
folder_name="Pan_troglodytes_PhyloP_MaskPanPanPanTro"

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2fasta /scratch/users/astarr97/PhyloP/hg38.447way.hal $phylop_species > $phylop_species.fasta
samtools faidx $phylop_species.fasta
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $phylop_species.fasta.fai | sort -n -r -k 3,3 > ${phylop_species}_contigs.bed

python make_scripts_mask_phylop.py ${phylop_species}_contigs.bed $to_mask $phylop_species $folder_name

cp /home/groups/hbfraser/astarr_scripts/AccelConv/driver.sh /scratch/users/astarr97/PhyloP/Feasibility/$folder_name
cd $folder_name
./driver.sh
cd ..
