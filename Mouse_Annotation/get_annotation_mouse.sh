#Downloaded the annotation from here for consistency's sake (V36)
#https://www.gencodegenes.org/mouse/
#Downloaded from ensembl mm10 website biomart human-mouse orthologs mm10_Ensembl102_hg38_orthologs.txt
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz
gunzip gencode.vM36.annotation.gtf.gz

#Get the cds and start codon exons, for now the start codons are unused
#Maybe need to change "gene_id" to what we actually want to name genes from the GTF.
python get_cds_start_reformat_mouse.py gencode.vM36.annotation.gtf gene_name

#Liftover to mm10
./liftOver gencode.vM36.annotation_CDS.bed mm39ToMm10.over.chain.gz gencode.vM36.annotation_CDS.mm10.bed gencode.vM36.annotation_CDS.mm10.err

sort -u gencode.vM36.annotation_CDS.mm10.bed > gencode.vM36.annotation_CDS.mm10.dedup.bed

#Need to merge CDS regions that partially overlap and are assigned to the same gene name
sort -k4,4 gencode.vM36.annotation_CDS.mm10.dedup.bed > gencode.vM36.annotation_CDS.mm10.sortgene.bed

#Split by gene
mkdir Dedup_CDS
python pergene_cds.py gencode.vM36.annotation_CDS.mm10.sortgene.bed

cd Dedup_CDS

#Merge bedfiles for each gene separately
for file in *.bed;
do
    sort -k1,1 -k2,2n $file > ${file::-4}.sort.bed
    bedtools merge -c 4 -o distinct -i ${file::-4}.sort.bed > ${file::-4}.merged.bed
done

#Make bed file with all genes
cat *.merged.bed > ../gencode.vM36.annotation_CDS.mm10.fulldedup.bed

cd ..

#Sort and get ready to intersect
sort -k1,1 -k2,2n gencode.vM36.annotation_CDS.mm10.fulldedup.bed > gencode.vM36.annotation_CDS.mm10.sort.bed

python make_human_gene_name.py