#human.gtf from one stored on OAK
#python get_cds_start_reformat_mouse.py human.gtf gene_name
#sort -u human_CDS.bed > human_CDS.uniq.bed
#sort -k4,4 human_CDS.uniq.bed > human_CDS.sortgene.bed

mkdir Dedup_CDS
python pergene_cds.py human_CDS.sortgene.bed

cd Dedup_CDS

#Merge bedfiles for each gene separately
for file in *.bed;
do
    sort -k1,1 -k2,2n $file > ${file::-4}.sort.bed
    bedtools merge -c 4 -o distinct -i ${file::-4}.sort.bed > ${file::-4}.merged.bed
done

#Make bed file with all genes
cat *.merged.bed > ../human_CDS.fulldedup.bed

cd ..

sort -k1,1 -k2,2n human_CDS.fulldedup.bed > /oak/stanford/groups/hbfraser/astarr/AccelConvDist/Annotations/human_CDS.sort.bed
