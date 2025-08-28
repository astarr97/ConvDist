o = open("to_rerun.txt")
out = open("rerun.sh", 'w')

prev_fold = 0 
for line in o:
    folder = line.split()[0]
    if folder != prev_fold and prev_fold:
        out.write("cd ../" + folder + "\n")
        out.write("cd All\nrm Run114.tsv\ncd ../run114\nrm *.out\nsbatch -p hbfraser run114.sh\ncd ..\n\n".replace("114", line.split()[1].replace("run", "").replace(".sh", "")))
    elif folder != prev_fold:
        out.write("cd " + folder + "\n")
        out.write("cd All\nrm Run114.tsv\ncd ../run114\nrm *.out\nsbatch -p hbfraser run114.sh\ncd ..\n\n".replace("114", line.split()[1].replace("run", "").replace(".sh", "")))
    else:
        out.write("cd All\nrm Run114.tsv\ncd ../run114\nrm *.out\nsbatch -p hbfraser run114.sh\ncd ..\n\n".replace("114", line.split()[1].replace("run", "").replace(".sh", "")))
    prev_fold = folder
o.close()
out.close()
