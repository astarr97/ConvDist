import sys

config = sys.argv[1]

species = []
o = open("hg38.447way.scientificNames.nh.txt")
for line in o:
    l = line.replace("(", "").replace(")", "").replace(" ", "").replace(",", "").split(":")[0]
    species.append(l)
o.close()

o = open(config)
c = 1
for line in o:
    if c == 1:
        assert(line=="Focal_species\tRelated_species\tOutgroup_species\tContigs_file\n")
    else:
        if len(line.split("\t")) != 4:
            print("Too few entries or separated by spaces instead of tabs in line " + str(c) + " of " + config)
            assert(False)
        else:
            l = line.split("\t")
            for entry in l[:-1]:
                entry_new = entry.split(";")
                for ent in entry_new:
                    e = ent.split(",")
                    for spec in e:
                        if spec not in species:
                            print(spec + " not in tree, likely spelling error")
                            assert(False)
            if "_contigs.bed" not in l[-1]:
                print("Must have contigs file as last entry")
                assert(False)
    c += 1
o.close()
