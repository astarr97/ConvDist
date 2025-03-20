orthos = pd.read_csv("mm10_Ensembl102_hg38_orthologs.txt", sep = "\t").dropna()
orthos["JustCaps"] = orthos["Gene name"].str.upper() == orthos["Human gene name"]
orthos = orthos.sort_values("Gene name")
d = {}
for index, row in orthos.iterrows():
    d[row["Gene name"]] = row["Human gene name"]

bed = pd.read_csv("gencode.vM36.annotation_CDS.mm10.sort.bed", sep = "\t", header = None)

c = 0
new_gene_name = []
for index, row in bed.iterrows():
    c += 1
    if c % 1000 == 0:
    if row[3] in d.keys():
        new_gene_name.append(d[row[3]])
    else:
        new_gene_name.append(row[3])

bed[3] = new_gene_name
bed.to_csv("gencode.vM36.annotation_CDS.mm10.hg38_name.sort.bed", sep = "\t", header = False, index = False)
