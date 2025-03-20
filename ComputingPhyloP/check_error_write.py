import os

out = open("move_error.sh", 'w')

prev_line = ""
prev_prev_line = ""
for file in os.listdir("./"):
    if ".out" in file:
        o = open(file)
        for line in o:
            if "ERROR" == line[0:5] and "ERRORED" not in prev_line:
                #Very important to comment out any print statements
                #print(prev_line, prev_prev_line)
                write_out = 0
                if prev_line.replace("\n", "") in os.listdir("./"):
                    out.write("mv " + prev_line.replace("\n", "") + " ERRORED_" + prev_line.replace("\n", "") + "\n\n")
                    out.write('file1=' + '"ERRORED_' + prev_line.replace("\n", "") + '"\n')
                    write_out = 1
                    contig_name = prev_line.replace("\n", "").split("_")[2]
                elif prev_prev_line.replace("\n", "") in os.listdir("./") and prev_line.replace("\n", "") not in os.listdir("./"):
                    out.write("mv " + prev_prev_line.replace("\n", "") + " ERRORED_" + prev_prev_line.replace("\n", "") + "\n\n")
                    out.write('file1=' + '"ERRORED_' + prev_prev_line.replace("\n", "") + '"\n')
                    write_out = 1
                    contig_name = prev_prev_line.replace("\n", "").split("_")[2]
                if write_out:
                    out.write("echo $file1\n")
                    out.write("/home/groups/hbfraser/Common_Software/phast/bin/maf_parse --split 1 --out-root ${file1::-4}_ERR $file1\n\n")
                    out.write("for file in ${file1::-4}_ERR*;\ndo\n\t")
                    out.write("/scratch/users/astarr97/Common_Software/phast/bin/phyloP --no-prune --chrom " + contig_name + " --msa-format MAF --method LRT --mode CONACC -d 6 --wig-scores -g /scratch/users/astarr97/astarr_scripts/AccelConvDist/fullTreeAnc239.100kb.mod $file > ${file::-4}.PhyloP.wig\n\t")
                    #out.write("python /scratch/users/astarr97/astarr_scripts/AccelConvDist/fix_wig.py ${file::-4}.PhyloP.wig " + contig_name +"\n\t")
                    out.write("/scratch/users/astarr97/astarr_scripts/AccelConvDist/bin/convert2bed --do-not-sort -i wig < ${file::-4}.PhyloP.wig > ${file::-4}.PhyloP.bed\n")
                    out.write("done\n\ncat ${file1::-4}_ERR*.PhyloP.bed > ${file1::-4}.PhyloP.bed\ncat ${file1::-4}_ERR*.PhyloP.bed > ${file1::-4}.PhyloP.bed\nrm ${file1::-4}_ERR*\n")
                    out.write("mv $file1 ../Errored\n\n")
            prev_prev_line = prev_line
            prev_line = line
        o.close()
out.close()
