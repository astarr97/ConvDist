#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB

#MUST BE IN BASE ENV TO HAVE RIGHT VERSION OF SCIPY FOR P-VALUE COMBINATION
#conda activate base

file_list="Orcinus_orca_Aquatic.FiltConv.PhyloP.NearestGene.SpecSup.bed.gz,Odobenus_rosmarus_Aquatic.FiltConv.Orcinus_orca.PhyloP.NearestGene.SpecSup.bed.gz"
folder_name="Aquatic"
spec_sup="100"
rem_single="0"
ref_species="Orcinus_orca"
python run_fwc_pwc.py $file_list $folder_name $spec_sup $rem_single $ref_species

file_list="Orcinus_orca_Aquatic.FiltConv.Rel.PhyloP.NearestGene.SpecSup.bed.gz,Odobenus_rosmarus_Aquatic.FiltConv.Rel.Orcinus_orca.PhyloP.NearestGene.SpecSup.bed.gz"
folder_name="Aquatic_Rel"
spec_sup="100"
rem_single="0"
ref_species="Orcinus_orca"
python run_fwc_pwc.py $file_list $folder_name $spec_sup $rem_single $ref_species

file_list="Orcinus_orca_Aquatic.FiltConv.PhyloP.NearestGene.SpecSup.bed.gz,Odobenus_rosmarus_Aquatic.FiltConv.Orcinus_orca.PhyloP.NearestGene.SpecSup.bed.gz"
folder_name="Aquatic_NoSingle"
spec_sup="100"
rem_single="1"
ref_species="Orcinus_orca"
python run_fwc_pwc.py $file_list $folder_name $spec_sup $rem_single $ref_species

file_list="Orcinus_orca_Aquatic.FiltConv.Rel.PhyloP.NearestGene.SpecSup.bed.gz,Odobenus_rosmarus_Aquatic.FiltConv.Rel.Orcinus_orca.PhyloP.NearestGene.SpecSup.bed.gz"
folder_name="Aquatic_Rel_NoSingle"
spec_sup="100"
rem_single="1"
ref_species="Orcinus_orca"
python run_fwc_pwc.py $file_list $folder_name $spec_sup $rem_single $ref_species
