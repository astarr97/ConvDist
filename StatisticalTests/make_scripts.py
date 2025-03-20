for i in range(1000):
    if i % 20 == 0:
        if i:
            out.close()
        out = open("run" + str(i//20) + ".sh", 'w')
        out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hns,hbfraser\n#SBATCH --mem=32GB\n\n")
        out.write('file_list="Orcinus_orca.Odobenus_rosmarus.Aquatic.ClosestVar.bed.gz,Orcinus_orca.Enhydra_lutris.Aquatic.ClosestVar.bed.gz,Enhydra_lutris.Odobenus_rosmarus.Aquatic.ClosestVar.bed.gz,Orcinus_orca.Trichechus_manatus.Aquatic.ClosestVar.bed.gz,Odobenus_rosmarus.Trichechus_manatus.Aquatic.ClosestVar.bed.gz,Enhydra_lutris.Trichechus_manatus.Aquatic.ClosestVar.bed.gz"\n\n')
        out.write('folder_name="Aquatic"\nprefix="Aquatic"\nspec_sup="100"\ndist_cut="1000"\nspecies_dif_cut="10"\npermute="1"\n\n')
        
    out.write('seed="' + str(i) + '"\n')
    out.write("python run_nearest_site.py $file_list $folder_name $prefix $spec_sup $dist_cut $species_dif_cut $permute $seed\n")
        
out.close()
