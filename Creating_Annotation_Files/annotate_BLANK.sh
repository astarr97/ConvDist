#For Nycticebus_pygmaeus we used the following modifications
#wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates/Otolemur_garnettii__small-eared_galago__otoGar3/geneAnnotation.gtf.gz
#gunzip geneAnnotation.gtf.gz
#python get_cds_start_reformat_TOGA.py geneAnnotation.gtf gene_id
#python add_point1.py geneAnnotation_CDS.bed
#/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal Otolemur_garnettii geneAnnotation_CDS.point1.bed Nycticebus_pygmaeus geneAnnotation_CDS.Nycticebus_pygmaeus.bed
#sort -u geneAnnotation_CDS.Nycticebus_pygmaeus.bed > geneAnnotation_CDS.Nycticebus_pygmaeus.uniq.bed
#mv geneAnnotation_CDS.Nycticebus_pygmaeus.uniq.bed geneAnnotation_CDS.dedup.bed

#For panda we are looking at https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002007445.2/ as the genome for the TOGA annotation
#But a different version of the genome was used for the hal file and we don't have liftover files between the two
#So instead we are going to use Ursus_maritimus which at least matches
#wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Carnivora/Ursus_maritimus__polar_bear__ursMar1/geneAnnotation.gtf.gz
#gunzip geneAnnotation.gtf.gz
#python get_cds_start_reformat_TOGA.py geneAnnotation.gtf gene_id
#python add_point1.py geneAnnotation_CDS.bed
#/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal Ursus_maritimus geneAnnotation_CDS.point1.bed Ailuropoda_melanoleuca geneAnnotation_CDS.Ailuropoda_melanoleuca.bed
#sort -u geneAnnotation_CDS.Ailuropoda_melanoleuca.bed > geneAnnotation_CDS.dedup.bed

#For Cricetulus_griseus, this is the version of the genome that is in the hal file https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_900186095.1/
#Which is not what is included as part of the TOGA output, that is the more recent version of the genome
#There is not suitable liftover unfortunately...
#So we use the Mesocricetus_auratus annotation and halLiftover
#Notably, the hibernation comparison does not include humans (only very diverged primates who debatably hibernate)
#Therefore, our strategy is to use the mouse TOGA
wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/mouse_mm10_reference/Rodentia/Mesocricetus_auratus__golden_hamster__mesAur1/geneAnnotation.gtf.gz
gunzip geneAnnotation.gtf.gz
python get_cds_start_reformat_TOGA.py geneAnnotation.gtf gene_id
python add_point1.py geneAnnotation_CDS.bed
/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal Mesocricetus_auratus geneAnnotation_CDS.point1.bed Cricetulus_griseus geneAnnotation_CDS.Cricetulus_griseus.bed
sort -u geneAnnotation_CDS.Cricetulus_griseus.bed > geneAnnotation_CDS.dedup.bed

#No modifications needed for CanFam4 or Echinops_telfairi or Dasypus novemcinctus
#wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Carnivora/Canis_lupus_familiaris__dog__canFam4/geneAnnotation.gtf.gz
#wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Afrotheria/Echinops_telfairi__small_Madagascar_hedgehog__echTel2/geneAnnotation.gtf.gz
#wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Xenarthra/Dasypus_novemcinctus__nine-banded_armadillo__dasNov3/geneAnnotation.gtf.gz

#Modifications made for Felis_catus_fca126, general strategy was to lift the annotation from FelCat9 to FelCat8, rename chromosome to match the hal file, then lift from FelCat8 (which is Felis_catus in the hal) to Felis_catus_fca126
#wget --no-check-certificate https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Carnivora/Felis_catus__domestic_cat__felCat9/geneAnnotation.gtf.gz
#gunzip geneAnnotation.gtf.gz
#python get_cds_start_reformat_TOGA.py geneAnnotation.gtf gene_id
#./liftOver -minMatch=0.95 geneAnnotation_CDS.bed felCat9ToFelCat8.over.chain.gz geneAnnotation_CDS.FelCat8.bed geneAnnotation_CDS.FelCat8.err
#Get the chromosome information for FelCat8

#wget https://hgdownload.soe.ucsc.edu/goldenPath/felCat8/database/chromAlias.txt.gz
#gunzip chromAlias.txt.gz
#python fix_cat_ChromNames.py
#Manually deleted the mitochondrial genes on chrom U20753.1
#/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/halLiftover --bedType 3 /scratch/users/astarr97/PhyloP/hg38.447way.hal Felis_catus geneAnnotation_CDS.FelCat8.renamed.bed Felis_catus_fca126 geneAnnotation_CDS.Felis_catus_fca126.bed
#sort -u geneAnnotation_CDS.Felis_catus_fca126.bed > geneAnnotation_CDS.dedup.bed

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
cat *.merged.bed > ../geneAnnotation_CDS_Felis_catus_fca126.fulldedup.bed

cd ..

#Sort and get ready to intersect
sort -k1,1 -k2,2n geneAnnotation_CDS_Felis_catus_fca126.fulldedup.bed > geneAnnotation_CDS_Felis_catus_fca126.sort.bed


