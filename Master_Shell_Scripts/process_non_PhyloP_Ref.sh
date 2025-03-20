#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -p hns,hbfraser
#SBATCH --mem=16GB

### IMPORTANT NOTE ###
#This will fail bedtools intersect for some species with annoyingly named scaffold-level genomes (despite me setting LC_COLLATE appropriately)
#The solution is to uncomment the cut -f 1,3 ../../${species}_contigs.bed > ${species}_contigs.txt line and the following line (around line 70)
#You then need to add the argument -g ${species}_contigs.sort.txt to the bedtools intersect call and that should solve it!

LC_COLLATE=C

#Define prefix and cat everything together
prefix="Odobenus_rosmarus_Aquatic"
species="Odobenus_rosmarus"
ref_species="Orcinus_orca"
phylop_path="/scratch/users/astarr97/PhyloP/Feasibility/Orcinus_orca_PhyloP_MaskAquatic/All/"
annotation_path="/oak/stanford/groups/hbfraser/astarr/AccelConv/Matriarchal/Orcinus_orca/"
scripts_path="/home/groups/hbfraser/astarr_scripts/AccelConvDist/"
output_dir="/scratch/users/astarr97/PhyloP/Feasibility/Aquatic_Convergence"
backup_path="/oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/Output_Files/"

mkdir /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}

cat *.tsv > ${prefix}_AllVariants.tsv
cat *.bed > ${prefix}_AllVariants_ForPhyloP.bed

#Sort and backup
sort -u ${prefix}_AllVariants_ForPhyloP.bed > ${prefix}_AllVariants_ForPhyloP.uniq.bed
sort -k1,1 -k2,2n ${prefix}_AllVariants_ForPhyloP.uniq.bed > ${prefix}_AllVariants_ForPhyloP.sort.bed
gzip -c ${prefix}_AllVariants_ForPhyloP.sort.bed > ${prefix}_AllVariants_ForPhyloP.sort.bed.gz
cp ${prefix}_AllVariants_ForPhyloP.sort.bed.gz /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}

sort -u ${prefix}_AllVariants.tsv > ${prefix}_AllVariants.dedup.tsv
sort -k1,1 -k2,2n ${prefix}_AllVariants.dedup.tsv > ${prefix}_AllVariants.sort.tsv
gzip ${prefix}_AllVariants.sort.tsv
cp ${prefix}_AllVariants.sort.tsv.gz /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}

mkdir Run_Files
mv *_ALL.tsv Run_Files
mv *_ALL_ForPhyloP.bed Run_Files

cp ${scripts_path}filter_for_conv.py ./
cp ${scripts_path}filter_for_poly.py ./
cp ${scripts_path}reformat_after_intersect.py ./
cp ${scripts_path}remove_dot.py ./

#Identify fixed and polymorphic sites within the clade (there may not be any polymorphic ones)
python filter_for_poly.py ${prefix}_AllVariants_ForPhyloP.sort.bed
python filter_for_conv.py ${prefix}_AllVariants_ForPhyloP.sort.bed

cp ${scripts_path}add_pos.py ./
cp ${scripts_path}get_bad_dups.py ./
cp ${scripts_path}filter_match_pos_fsamepos.py ./

#Start the reciprocal liftover process to only get one-to-one sites
python add_pos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.bed
python add_pos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.bed

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.pos.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.FOC.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.FOC.bed ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.bed

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.pos.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.FOC.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.FOC.bed ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.bed

#This version does not remove sites that do not lift over one-to-one
#python filter_match_pos.py All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.REL.FOC.bed All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.lifted.bed
#sort -k1,1 -k2,2n All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.lifted.bed > All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.sort.bed

#cut -f 1,3 ../../${species}_contigs.bed > ${species}_contigs.txt
#sort -k1,1 ${species}_contigs.txt > ${species}_contigs.sort.txt

#This version does remove those sites
sort -k7,7 ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.sortforbad.bed
python get_bad_dups.py ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.sortforbad.bed ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.bad.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.bad.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.bad.sort.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.sort.bed
bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.sort.bed -b ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.bad.sort.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.ToRemove.bed
python filter_match_pos_fsamepos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.BACK.ToRemove.bed ${prefix}_AllVariants_ForPhyloP.FiltConv.ToLift.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.FiltConv.ToLift.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.FiltConv.${ref_species}.bed
sort -k1,1 -k2,2n ${prefix}_AllVariants_ForPhyloP.FiltConv.${ref_species}.bed > ${prefix}_AllVariants_ForPhyloP.FiltConv.${ref_species}.sort.bed

#Repeat filtering process for Rel species variants
sort -k7,7 ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.sortforbad.bed
python get_bad_dups.py ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.sortforbad.bed ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.bad.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.bad.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.bad.sort.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.sort.bed
bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.sort.bed -b ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.bad.sort.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.ToRemove.bed
python filter_match_pos_fsamepos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.BACK.ToRemove.bed ${prefix}_AllVariants_ForPhyloP.FiltConv.ToLift.Rel.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.FiltConv.ToLift.Rel.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.FiltConv.Rel.${ref_species}.bed
sort -k1,1 -k2,2n ${prefix}_AllVariants_ForPhyloP.FiltConv.Rel.${ref_species}.bed > ${prefix}_AllVariants_ForPhyloP.FiltConv.Rel.${ref_species}.sort.bed


#Do for the polymorphic sites (if they exist)
#Start the reciprocal liftover process to only get one-to-one sites
python add_pos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.bed
python add_pos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.bed

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.pos.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.FOC.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.FOC.bed ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.bed

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.pos.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.FOC.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${ref_species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.FOC.bed ${species} ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.bed

#This version does not remove sites that do not lift over one-to-one
#python filter_match_pos.py All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.REL.FOC.bed All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.lifted.bed
#sort -k1,1 -k2,2n All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.lifted.bed > All_${ref_species}_PhyloP_Dedup_FiltCut$PhyloP_Cutoff.sort.bed

#This version does remove those sites
sort -k7,7 ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.sortforbad.bed
python get_bad_dups.py ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.sortforbad.bed ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.bad.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.bad.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.bad.sort.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.sort.bed
bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.sort.bed -b ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.bad.sort.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.ToRemove.bed
python filter_match_pos_fsamepos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.BACK.ToRemove.bed ${prefix}_AllVariants_ForPhyloP.FiltPoly.ToLift.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.FiltPoly.ToLift.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.FiltPoly.${ref_species}.bed
sort -k1,1 -k2,2n ${prefix}_AllVariants_ForPhyloP.FiltPoly.${ref_species}.bed > ${prefix}_AllVariants_ForPhyloP.FiltPoly.${ref_species}.sort.bed

#This version does remove those sites
sort -k7,7 ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.sortforbad.bed
python get_bad_dups.py ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.sortforbad.bed ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.bad.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.bad.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.bad.sort.bed
sort -u ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.bed | sort -k1,1 -k2,2n > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.sort.bed
bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.sort.bed -b ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.bad.sort.bed > ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.ToRemove.bed
python filter_match_pos_fsamepos.py ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.BACK.ToRemove.bed ${prefix}_AllVariants_ForPhyloP.FiltPoly.ToLift.Rel.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}_AllVariants_ForPhyloP.FiltPoly.ToLift.Rel.bed ${ref_species} ${prefix}_AllVariants_ForPhyloP.FiltPoly.Rel.${ref_species}.bed
sort -k1,1 -k2,2n ${prefix}_AllVariants_ForPhyloP.FiltPoly.Rel.${ref_species}.bed > ${prefix}_AllVariants_ForPhyloP.FiltPoly.Rel.${ref_species}.sort.bed

cp ${scripts_path}reformat_after_intersect.py ./
cp ${scripts_path}remove_dot.py ./

#Intersect with PhyloP, nearest gene, and species support
bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.FiltConv.${ref_species}.sort.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquatic.sort.NoID.bed > ${prefix}.FiltConv.${ref_species}.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltConv.${ref_species}.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_${ref_species}_CDS.uniq.sort.bed > ${prefix}.FiltConv.${ref_species}.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltConv.${ref_species}.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltConv.${ref_species}.PhyloP.NearestGene.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquatic.SpecSup.sort.bed > ${prefix}.FiltConv.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed

bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.FiltPoly.${ref_species}.sort.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquatic.sort.NoID.bed > ${prefix}.FiltPoly.${ref_species}.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltPoly.${ref_species}.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_${ref_species}_CDS.uniq.sort.bed > ${prefix}.FiltPoly.${ref_species}.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltPoly.${ref_species}.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltPoly.${ref_species}.PhyloP.NearestGene.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquatic.SpecSup.sort.bed > ${prefix}.FiltPoly.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed

#Repeat for rel species
#Need to change PhyloP path and add Flippedto the names of the PhyloP and SpecSup files
phylop_path="/scratch/users/astarr97/PhyloP/Feasibility/Orcinus_orca_PhyloP_MaskAquaticFlipped/All/"

bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.FiltConv.Rel.${ref_species}.sort.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquaticFlipped.sort.NoID.bed > ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_${ref_species}_CDS.uniq.sort.bed > ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.NearestGene.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquaticFlipped.SpecSup.sort.bed > ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed

bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.FiltPoly.Rel.${ref_species}.sort.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquaticFlipped.sort.NoID.bed > ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_${ref_species}_CDS.uniq.sort.bed > ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.NearestGene.bed -b ${phylop_path}All_${ref_species}_PhyloP_MaskAquaticFlipped.SpecSup.sort.bed > ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed

#Fix it up so that the extraneous columns are deleted
python reformat_after_intersect.py ${prefix}.FiltConv.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed
python reformat_after_intersect.py ${prefix}.FiltPoly.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed

#Copy to the directory where the analysis will continue
cp ${prefix}.FiltConv.${ref_species}.PhyloP.NearestGene.SpecSup.bed $output_dir
cp ${prefix}.FiltPoly.${ref_species}.PhyloP.NearestGene.SpecSup.bed $output_dir

#Repeat for rel species
python reformat_after_intersect.py ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed
python reformat_after_intersect.py ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.NearestGene.SpecSup.ToFix.bed

cp ${prefix}.FiltConv.Rel.${ref_species}.PhyloP.NearestGene.SpecSup.bed $output_dir
cp ${prefix}.FiltPoly.Rel.${ref_species}.PhyloP.NearestGene.SpecSup.bed $output_dir

cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/get_singeltons.py ./

#Finally, extract the singelton information and back it all up
python get_singeltons.py ${prefix}_AllVariants_ForPhyloP.sort.bed $prefix

/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}.Singelton.Focal.ToFix.bed ${ref_species} ${prefix}.Singelton.Focal.${ref_species}.ToFix.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}.Singelton.Rel.ToFix.bed ${ref_species} ${prefix}.Singelton.Rel.${ref_species}.ToFix.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}.Singelton.Out.ToFix.bed ${ref_species} ${prefix}.Singelton.Out.${ref_species}.ToFix.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}.MultiSingelton.Focal.ToFix.bed ${ref_species} ${prefix}.MultiSingelton.Focal.${ref_species}.ToFix.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}.MultiSingelton.Rel.ToFix.bed ${ref_species} ${prefix}.MultiSingelton.Rel.${ref_species}.ToFix.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal ${species} ${prefix}.MultiSingelton.Out.ToFix.bed ${ref_species} ${prefix}.MultiSingelton.Out.${ref_species}.ToFix.bed

mkdir /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot

gzip -c ${prefix}.Singelton.Focal.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Focal.ToFix.bed.gz
gzip -c ${prefix}.Singelton.Rel.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Rel.ToFix.bed.gz
gzip -c ${prefix}.Singelton.Out.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Out.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Focal.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Focal.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Rel.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Rel.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Out.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Out.ToFix.bed.gz

gzip -c ${prefix}.Singelton.Focal.${ref_species}.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Focal.${ref_species}.ToFix.bed.gz
gzip -c ${prefix}.Singelton.Rel.${ref_species}.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Rel.${ref_species}.ToFix.bed.gz
gzip -c ${prefix}.Singelton.Out.${ref_species}.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Out.${ref_species}.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Focal.${ref_species}.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Focal.${ref_species}.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Rel.${ref_species}.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Rel.${ref_species}.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Out.${ref_species}.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Out.${ref_species}.ToFix.bed.gz

for file in *Singelton*ToFix.bed;
do
    python remove_dot.py $file
done

gzip -c ${prefix}.Singelton.Focal.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.Singelton.Focal.bed.gz
gzip -c ${prefix}.Singelton.Rel.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.Singelton.Rel.bed.gz
gzip -c ${prefix}.Singelton.Out.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.Singelton.Out.bed.gz
gzip -c ${prefix}.MultiSingelton.Focal.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.MultiSingelton.Focal.bed.gz
gzip -c ${prefix}.MultiSingelton.Rel.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.MultiSingelton.Rel.bed.gz
gzip -c ${prefix}.MultiSingelton.Out.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.MultiSingelton.Out.bed.gz

gzip -c ${prefix}.Singelton.Focal.${ref_species}.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.Singelton.Focal.${ref_species}.bed.gz
gzip -c ${prefix}.Singelton.Rel.${ref_species}.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.Singelton.Rel.${ref_species}.bed.gz
gzip -c ${prefix}.Singelton.Out.${ref_species}.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.Singelton.Out.${ref_species}.bed.gz
gzip -c ${prefix}.MultiSingelton.Focal.${ref_species}.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.MultiSingelton.Focal.${ref_species}.bed.gz
gzip -c ${prefix}.MultiSingelton.Rel.${ref_species}.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.MultiSingelton.Rel.${ref_species}.bed.gz
gzip -c ${prefix}.MultiSingelton.Out.${ref_species}.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/${prefix}.MultiSingelton.Out.${ref_species}.bed.gz

