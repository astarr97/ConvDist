import sys

#First argument is the bed file from bedtools getfasta, second is the reference genome species name
fasta_bed = sys.argv[1]
species = sys.argv[2]

assert(fasta_bed.endswith(".fastabed"))

#The only difference between the two output files should be that the central base in the sequence is different, there is a lot of redundancy
o = open(fasta_bed)
out1 = open(fasta_bed.replace(".Comped.2114.fastabed", "." + species + ".Final.fastabed"), 'w')
out2 = open(fasta_bed.replace(".Comped.2114.fastabed", "." + species + ".ToCompare.fastabed"), 'w')

rev = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}

c = 0
t = 0
#Iterate through each line
for line in o:
    t += 1
    l = line.replace("\n", "").split("\t")
    info = l[3].split(";")
    
    #Position of the base to change
    position = 1056
    seq = list(l[4].upper())
    
    #Pull the ancestral base to make sure we are doing things on the right strand
    if "Focal" in fasta_bed:
        assert_base = info[5]
    elif "Closest" in fasta_bed:
        assert_base = info[3]
    
    #Check if we need to take the complement of it by comparing the reference base from the TSV to the reference base in the sequence from getfasta
    if rev[seq[position]] == assert_base:
        new_info = [info[0]]
        for i in info[1:]:
            new_info.append(rev[i])
        new_info.append("Rev")
    elif seq[position] == assert_base:
        new_info = info
        new_info.append("NotRev")
    #If it is not the correct base or the complement, then print error and exit
    else:
        print("Ancestral base mismatch", seq[1056], assert_base)
        assert(False)

    #Assume that the first species in the comparison is listed first, and second is listed second separated by a "."
    if fasta_bed.split(".")[0] == species and "Focal" in fasta_bed:
        new_base_foc = new_info[1]
        new_base_rel = new_info[3]
    elif fasta_bed.split(".")[1] == species and "Focal" in fasta_bed:
        new_base_foc = new_info[2]
        new_base_rel = new_info[4]
    elif fasta_bed.split(".")[0] == species and "Closest" in fasta_bed:
        new_base_foc = new_info[1]
        new_base_rel = new_info[2]
    elif fasta_bed.split(".")[1] == species and "Closest" in fasta_bed:
        new_base_foc = new_info[1]
        new_base_rel = new_info[2]
    
    #Check to make sure that the rel and focal bases are not equal, they should never be equal
    try:
        assert(new_base_foc != new_base_rel)
    except:
        #print(new_base_foc, new_base_rel, "Rel and focal are equal suggesting error")
        #print(line)
        #assert(False)
        c += 1
    
    #Make changes and write out
    seq[1056] = new_base_foc
    out1.write("\t".join(l[0:3] + [";".join(new_info), "".join(seq)]) + "\n")
    
    seq[1056] = new_base_rel
    out2.write("\t".join(l[0:3] + [";".join(new_info), "".join(seq)]) + "\n")

if c:
    print(c, t, c/t, "Some sites had rel equal to foc, if it is a non-negligible number there was likely an error", fasta_bed)
out1.close()
out2.close()
o.close()
