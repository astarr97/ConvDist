#Define tsv_file of interest
tsv_file="Helogale_parvula_Matriarchal_AllVariants.sort.tsv.gz"
header_file="ForCooperativeBreeding_Helogale_parvula_header.txt"


#Define focal, related, and outgroup species as well as species to add
focal="Helogale_parvula"
rel="Cryptoprocta_ferox"
outgroup="Paradoxurus_hermaphroditus"
add_focal="Suricata_suricatta"
add_rel="0"
add_out="0"

lift="0"

prefix="Helogale_parvula_ForCooperativeBreeding_"

#Copy needed script and halLiftover (copy over to prevent CooperativeBreedingdown due to I/O from scratch, home, and oak)
cp /scratch/users/astarr97/astarr_scripts/AccelConvDist/filter_and_tag_variants.py ./
cp /home/groups/hbfraser/Common_Software/cactus-bin-v2.6.13/bin/halLiftover ./

#Uncompress the file
zcat $tsv_file > ${tsv_file::-3}

#Get the species in the file
grep ref ${tsv_file::-3} > $header_file

#Split it into 200 megabyte chunks
split --line-bytes 200000000 ${tsv_file::-3} $prefix

for file in ${prefix}*;
do
    cat $header_file $file > $file.tsv
    rm $file
done

python make_scripts_pull_species.py $prefix $lift $focal $rel $outgroup $add_focal $add_rel $add_out

for file in run*.sh;
do
    sbatch -p hns,hbfraser $file
done

#Once that completes, cat everything together and sort just as if we had pulled the variants the normal way!
