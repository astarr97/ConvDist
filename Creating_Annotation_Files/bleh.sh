cd Dedup_CDS

#Merge bedfiles for each gene separately
for file in *.bed;
do
    sort -k1,1 -k2,2n $file > ${file::-4}.sort.bed
    bedtools merge -c 4 -o distinct -i ${file::-4}.sort.bed > ${file::-4}.merged.bed
done

#Make bed file with all genes
cat *.merged.bed > ../geneAnnotation_CDS_Dasypus_novemcinctus.fulldedup.bed

cd ..

#Sort and get ready to intersect
sort -k1,1 -k2,2n geneAnnotation_CDS_Dasypus_novemcinctus.fulldedup.bed > geneAnnotation_CDS_Dasypus_novemcinctus.sort.bed


