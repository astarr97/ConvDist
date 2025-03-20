#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB

LC_COLLATE=C

#Prefix to name the PhyloP
prefix="All_Orcinus_orca_PhyloP_MaskAquatic"
#Species to which the PhyloP is referenced
focal_species="Orcinus_orca"
#Path to back things up to
backup_path="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/"

mkdir ${backup_path}/PhyloP

cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/remove_id.py ./

#Copy needed scripts over
cp ../../${focal_species}_contigs.bed ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/remove_id.py ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/bedGraphToBigWig ./
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/spec_sup_make_bed.py ./

#Cat together, remove id column, convert to bigWig to back up
cat *.PhyloP.bed > $prefix.bed
sort -u $prefix.bed > $prefix.uniq.bed
sort -k1,1 -k2,2n $prefix.uniq.bed > $prefix.sort.bed

python remove_id.py $prefix.sort.bed
cut -f 1,3 ${focal_species}_contigs.bed > ${focal_species}_contigs.txt
./bedGraphToBigWig $prefix.sort.NoID.bed ${focal_species}_contigs.txt $prefix.bw
mv $prefix.bw ${backup_path}PhyloP

#Take another stab at backing up because the bigwig thing errors sometimes and I'm not sure it worked
gzip -c $prefix.sort.NoID.bed > ${backup_path}PhyloP/$prefix.sort.NoID.bed.gz

#Cat together species support, convert to bed, gzip, and backup
cat *.SpecSup.txt > $prefix.SpecSup.txt
python spec_sup_make_bed.py $prefix.SpecSup.txt
sort -k1,1 -k2,2n $prefix.SpecSup.bed > $prefix.SpecSup.sort.bed
gzip -c $prefix.SpecSup.sort.bed > $prefix.SpecSup.sort.bed.gz
mv $prefix.SpecSup.sort.bed.gz ${backup_path}PhyloP
