import os

#Iterate through the new TSV folders, check that there are enough files and that none of the logs have error messages
for folder in os.listdir():
    if "." not in folder:
        folders = os.listdir(folder)
        files = os.listdir(folder + "/All")
        print(len(folders), len(files), folder)
        if len(folders) -2 != len(files):
            print("Not enough files, " + folder)
        
        for fold in folders:
            if "." not in fold:
                for f in os.listdir(folder + "/" + fold):
                    if ".out" in f:
                        o = open(folder + "/" + fold + "/" + f)
                        for line in o:
                            if "rror" in line or "RROR" in line or "arning" in line or "egmentation" in line or "EGMENTATION" in line or "ARNING" in line or "exception" in line:
                                print(folder, fold, f, "Contains error")
