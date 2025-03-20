#Downloaded the annotation from here for consistency's sake (V36)
#https://www.gencodegenes.org/mouse/
#Downloaded from ensembl mm10 website biomart human-mouse orthologs mm10_Ensembl102_hg38_orthologs.txt
#Downloaded annotation from: https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Rodentia/Heterocephalus_glaber__naked_mole-rat__hetGla2/
gunzip geneAnnotation.gtf.gz

#Get the cds and start codon exons, for now the start codons are unused
#Maybe need to change "gene_id" to what we actually want to name genes from the GTF.
python get_cds_start_reformat_TOGA.py geneAnnotation.gtf gene_id


sort -u geneAnnotation_CDS.bed > geneAnnotation_CDS.dedup.bed

#Need to merge CDS regions that partially overlap and are assigned to the same gene name
sort -k4,4 geneAnnotation_CDS.dedup.bed > geneAnnotation_CDS.sortgene.bed

#Split by gene
mkdir Dedup_CDS
python pergene_cds.py geneAnnotation_CDS.sortgene.bed

cd Dedup_CDS

#Merge bedfiles for each gene separately
for file in *.bed;
do
    sort -k1,1 -k2,2n $file > ${file::-4}.sort.bed
    bedtools merge -c 4 -o distinct -i ${file::-4}.sort.bed > ${file::-4}.merged.bed
done

#Make bed file with all genes
cat *.merged.bed > ../geneAnnotation_CDS_Hetercephalus_glaber.fulldedup.bed

cd ..

#Sort and get ready to intersect
sort -k1,1 -k2,2n geneAnnotation_CDS_Hetercephalus_glaber.fulldedup.bed > geneAnnotation_CDS_Hetercephalus_glaber.sort.bed
