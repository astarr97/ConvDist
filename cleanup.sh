folder="Orcinus_orca_PhyloP_MaskAquatic_New"
cd $folder

for fold in run*;
do
    cd $fold
    rm *.maf
    cd ..
done
