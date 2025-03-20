#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH -p hns,hbfraser
#SBATCH --mem=32GB

LC_COLLATE=C

#Define prefix and cat everything together
prefix="Orcinus_orca_Aquatic"
species="Orcinus_orca"
phylop_path="/scratch/users/astarr97/PhyloP/Feasibility/Orcinus_orca_PhyloP_MaskAquatic/All/"
annotation_path="/oak/stanford/groups/hbfraser/astarr/AccelConv/Matriarchal/Orcinus_orca/"
scripts_path="/home/groups/hbfraser/astarr_scripts/AccelConvDist/"
output_dir="/scratch/users/astarr97/PhyloP/Feasibility/Aquatic_Convergence_New"

mkdir /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}

#Cat everything together
cat *.tsv > ${prefix}_AllVariants.tsv
cat *.bed > ${prefix}_AllVariants_ForPhyloP.bed

#Sort and backup
sort -u ${prefix}_AllVariants_ForPhyloP.bed > ${prefix}_AllVariants_ForPhyloP.uniq.bed
sort -k1,1 -k2,2n ${prefix}_AllVariants_ForPhyloP.uniq.bed > ${prefix}_AllVariants_ForPhyloP.sort.bed
gzip -c ${prefix}_AllVariants_ForPhyloP.sort.bed > ${prefix}_AllVariants_ForPhyloP.sort.bed.gz
cp ${prefix}_AllVariants_ForPhyloP.sort.bed.gz /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/Orcinus_orca

#sort -u ${prefix}_AllVariants.tsv > ${prefix}_AllVariants.dedup.tsv
sort -k1,1 -k2,2n ${prefix}_AllVariants.dedup.tsv > ${prefix}_AllVariants.sort.tsv
gzip ${prefix}_AllVariants.sort.tsv
cp ${prefix}_AllVariants.sort.tsv.gz /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/Orcinus_orca

#mkdir Run_Files
mv *_ALL.tsv Run_Files
mv *_ALL_ForPhyloP.bed Run_Files

cp ${scripts_path}filter_for_conv.py ./
cp ${scripts_path}filter_for_poly.py ./
cp ${scripts_path}reformat_after_intersect.py ./
cp ${scripts_path}remove_dot.py ./

#Identify fixed (for conv) and polymorphic within the clade (for poly) sites
python filter_for_poly.py ${prefix}_AllVariants_ForPhyloP.sort.bed
python filter_for_conv.py ${prefix}_AllVariants_ForPhyloP.sort.bed

#Intersect with PhyloP, nearest gene, and species support
bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquatic.sort.NoID.bed > ${prefix}.FiltConv.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltConv.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_Orcinus_orca_CDS.uniq.sort.bed > ${prefix}.FiltConv.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltConv.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltConv.PhyloP.NearestGene.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquatic.SpecSup.sort.bed > ${prefix}.FiltConv.PhyloP.NearestGene.SpecSup.ToFix.bed

bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquatic.sort.NoID.bed > ${prefix}.FiltPoly.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltPoly.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_Orcinus_orca_CDS.uniq.sort.bed > ${prefix}.FiltPoly.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltPoly.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltPoly.PhyloP.NearestGene.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquatic.SpecSup.sort.bed > ${prefix}.FiltPoly.PhyloP.NearestGene.SpecSup.ToFix.bed

#Fix it up so that the extraneous columns are deleted
python reformat_after_intersect.py ${prefix}.FiltConv.PhyloP.NearestGene.SpecSup.ToFix.bed
python reformat_after_intersect.py ${prefix}.FiltPoly.PhyloP.NearestGene.SpecSup.ToFix.bed

#Copy to the directory where the analysis will continue
cp ${prefix}.FiltConv.PhyloP.NearestGene.SpecSup.bed $output_dir
cp ${prefix}.FiltPoly.PhyloP.NearestGene.SpecSup.bed $output_dir

#Repeat everything for rel
#Need to change phylop_path and add Flipped at the end of the PhyloP and SpecSup files
phylop_path="/scratch/users/astarr97/PhyloP/Feasibility/Orcinus_orca_PhyloP_MaskAquaticFlipped/All/"

bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltConv.Rel.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquaticFlipped.sort.NoID.bed > ${prefix}.FiltConv.Rel.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltConv.Rel.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_Orcinus_orca_CDS.uniq.sort.bed > ${prefix}.FiltConv.Rel.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltConv.Rel.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltConv.Rel.PhyloP.NearestGene.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquaticFlipped.SpecSup.sort.bed > ${prefix}.FiltConv.Rel.PhyloP.NearestGene.SpecSup.ToFix.bed

bedtools intersect -sorted -wao -a ${prefix}_AllVariants_ForPhyloP.sort.FiltPoly.Rel.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquaticFlipped.sort.NoID.bed > ${prefix}.FiltPoly.Rel.PhyloP.bed
bedtools closest -d -wao -a ${prefix}.FiltPoly.Rel.PhyloP.bed -b ${annotation_path}geneAnnotation_Tursiops_truncatus_Lifted_Orcinus_orca_CDS.uniq.sort.bed > ${prefix}.FiltPoly.Rel.PhyloP.NearestGene.ToFix.bed
python remove_dot.py ${prefix}.FiltPoly.Rel.PhyloP.NearestGene.ToFix.bed
bedtools intersect -sorted -wao -a ${prefix}.FiltPoly.Rel.PhyloP.NearestGene.bed -b ${phylop_path}All_Orcinus_orca_PhyloP_MaskAquaticFlipped.SpecSup.sort.bed > ${prefix}.FiltPoly.Rel.PhyloP.NearestGene.SpecSup.ToFix.bed

#Fix it up so that the extraneous columns are deleted
python reformat_after_intersect.py ${prefix}.FiltConv.Rel.PhyloP.NearestGene.SpecSup.ToFix.bed
python reformat_after_intersect.py ${prefix}.FiltPoly.Rel.PhyloP.NearestGene.SpecSup.ToFix.bed

#Copy to the directory where the analysis will continue
cp ${prefix}.FiltConv.Rel.PhyloP.NearestGene.SpecSup.bed $output_dir
cp ${prefix}.FiltPoly.Rel.PhyloP.NearestGene.SpecSup.bed $output_dir

#Get the singelton information and backup
cp /home/groups/hbfraser/astarr_scripts/AccelConvDist/get_singeltons.py ./

python get_singeltons.py ${prefix}_AllVariants_ForPhyloP.sort.bed $prefix

mkdir /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot

gzip -c ${prefix}.Singelton.Focal.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Focal.ToFix.bed.gz
gzip -c ${prefix}.Singelton.Rel.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Rel.ToFix.bed.gz
gzip -c ${prefix}.Singelton.Out.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.Singelton.Out.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Focal.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Focal.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Rel.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Rel.ToFix.bed.gz
gzip -c ${prefix}.MultiSingelton.Out.ToFix.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Aquatic/${species}/WithDot/${prefix}.MultiSingelton.Out.ToFix.bed.gz

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
