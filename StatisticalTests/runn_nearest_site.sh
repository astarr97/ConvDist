#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=32GB

#MUST BE IN BASE ENV TO HAVE RIGHT VERSION OF SCIPY FOR P-VALUE COMBINATION
#conda activate base

file_list="Orcinus_orca.Odobenus_rosmarus.Aquatic.ClosestVar.bed.gz,Orcinus_orca.Enhydra_lutris.Aquatic.ClosestVar.bed.gz,Enhydra_lutris.Odobenus_rosmarus.Aquatic.ClosestVar.bed.gz,Orcinus_orca.Trichechus_manatus.Aquatic.ClosestVar.bed.gz,Odobenus_rosmarus.Trichechus_manatus.Aquatic.ClosestVar.bed.gz,Enhydra_lutris.Trichechus_manatus.Aquatic.ClosestVar.bed.gz"
folder_name="Aquatic"
prefix="Aquatic"
spec_sup="100"
dist_cut="1000"
species_dif_cut="10"
permute="0"
seed="6"
#python run_nearest_site.py $file_list $folder_name $prefix $spec_sup $dist_cut $species_dif_cut $permute

permute="1"
#python run_nearest_site.py $file_list $folder_name $prefix $spec_sup $dist_cut $species_dif_cut $permute $seed

#Removing this cuz we can't do it for the other one Enhydra_lutris.Odobenus_rosmarus.Aquatic.ClosestVar.bed.gz,
file_list="Orcinus_orca.Odobenus_rosmarus.Aquatic.ClosestVar.bed.gz,Orcinus_orca.Enhydra_lutris.Aquatic.ClosestVar.bed.gz,Orcinus_orca.Trichechus_manatus.Aquatic.ClosestVar.bed.gz,Odobenus_rosmarus.Trichechus_manatus.Aquatic.ClosestVar.bed.gz,Enhydra_lutris.Trichechus_manatus.Aquatic.ClosestVar.bed.gz"
folder_name="Aquatic_NoEnhOdo"
prefix="Aquatic"
spec_sup="100"
dist_cut="1000"
species_dif_cut="10"
permute="0"
seed="6"
#python run_nearest_site.py $file_list $folder_name $prefix $spec_sup $dist_cut $species_dif_cut $permute

permute="1"
#python run_nearest_site.py $file_list $folder_name $prefix $spec_sup $dist_cut $species_dif_cut $permute $seed


#Important to have Rel in files names if they are Rel!
#NOTE we cannot run Enhydra_lutris.Odobenus_rosmarus.Aquatic.Rel.ClosestVar.bed.gz because it has Mustela_putorius as the focal species for Enhydra_lutris Rel and as an outgroup for Odobenus_rosmarus Rel!
file_list="Orcinus_orca.Odobenus_rosmarus.Aquatic.Rel.ClosestVar.bed.gz,Orcinus_orca.Enhydra_lutris.Aquatic.Rel.ClosestVar.bed.gz,Orcinus_orca.Trichechus_manatus.Aquatic.Rel.ClosestVar.bed.gz,Odobenus_rosmarus.Trichechus_manatus.Aquatic.Rel.ClosestVar.bed.gz,Enhydra_lutris.Trichechus_manatus.Aquatic.Rel.ClosestVar.bed.gz"

folder_name="Aquatic_Rel"
spec_sup="100"
dist_cut="1000"
prefix="Aquatic.Rel"
species_dif_cut="10"
permute="0"
seed="6"
python run_nearest_site.py $file_list $folder_name $prefix $spec_sup $dist_cut $species_dif_cut $permute

permute="1"
python run_nearest_site.py $file_list $folder_name $prefix $spec_sup $dist_cut $species_dif_cut $permute $seed
