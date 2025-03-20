import pandas as pd

gene_names = pd.read_csv("chromAlias.txt", sep = '\t', header = None)
gene_names = gene_names[gene_names[2].isin(["genbank", "genbank,ensembl"])]
gene_names = gene_names[[0, 1]].set_index(1)
gene_names.columns = ["NewChromName"]

bed_file = pd.read_csv("geneAnnotation_CDS.FelCat8.bed", header = None, sep = "\t").set_index(0)
new_bed = bed_file.join(gene_names)

new_bed = new_bed[["NewChromName", 1, 2, 3]]
new_bed.to_csv("geneAnnotation_CDS.FelCat8.renamed.bed", sep = "\t", index = False, header = False)
