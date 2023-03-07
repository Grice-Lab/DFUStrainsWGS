# Amy Campbell
# Feb 2023
# Make a shellscript for each patient with a CSV in input csvs folder for cladebreaker
# also make a helpful little text file of bsub commands to kick off these jobs :)
import os

PathToInputs="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/InputCSVs/"
configstring="/home/acampbe/DFU/data/WGS_2020/cladebreaker/penn_cluster.config"


listInputs = os.listdir(PathToInputs)
listInputs.remove(".DS_Store")
#len(listInputs)
jobcommands=[]

OldJobIDname = listInputs[0].replace("_input.csv", "") + "_job"

outputfolders=[]
listInputs = list(filter(lambda k: 'CC' in k, listInputs))
for l in range(len(listInputs)):
    inputcsvpath = "/home/acampbe/DFU/data/WGS_2020/cladebreaker/InputCSVs/" + str(listInputs[l])
    jobIDstring=listInputs[l].replace("_input.csv", "")
    outputpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_" + jobIDstring
    inputpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/InputCSVs/"  + listInputs[l]
    NumberID=jobIDstring.replace("CC","")
    BatchFileName="Run_Cladebreaker_CC" + str(NumberID) +".sh"
    BatchFilePath = "/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/BatchShells/" + BatchFileName
    jobIDname= jobIDstring+ "_job"
    if jobIDname==OldJobIDname:
        jobRunCommand = "bsub -e " + jobIDstring + ".e -o " + jobIDstring + ".o" + " -J \"" + jobIDname +"\" sh " +  BatchFileName
    else:
        jobRunCommand= "bsub -e " + jobIDstring + ".e -o " + jobIDstring + ".o " + "-w \"done(" +OldJobIDname +")\" " + "-J \"" + jobIDname +"\" sh " +  BatchFileName

    OldJobIDname=jobIDname
    with open(BatchFilePath, "w") as outputshell:
        outputshell.write("#!bin/bash\n")
        outputshell.write("# cladebreaker\n")
        outputshell.write("################\n")
        outputshell.write("source ~/mambaforge/bin/activate ~/mambaforge/envs/cladebreaker2\n")
        outputshell.write("\n")
        outputshell.write("mkdir -p " + outputpath + "\n")
        outputshell.write("\n")
        outputshell.write("cladebreaker \\\n")
        outputshell.write("--input "+inputpath+" \\\n")
        outputshell.write("--outdir " + outputpath + " \\\n")
        outputshell.write("--coverage 100 --db /home/acampbe/DownloadedDatabases/WhatsGNU_Sau_Ortholog/Sau_Ortholog_10350.pickle \\\n")
        outputshell.write("--o --topgenomes_count 5 --force -profile conda -with-conda true \\\n")
        outputshell.write("-c /home/acampbe/DFU/data/WGS_2020/cladebreaker/penn_cluster.config")
        outputshell.write("\n")
        outputshell.write("\n")
        outputshell.write("# raxML\n")
        outputshell.write("################\n")
        outputshell.write("source ~/mambaforge/bin/activate ~/mambaforge/envs/TreeEnv2\n")
        outputshell.write("# Estimate tree based on roary output\n")
        outputshell.write("raxmlHPC -m GTRGAMMA -p 19104 \\\n")
        outputshell.write("-s " + outputpath+ "/roary_alignment/results/core_gene_alignment.aln" + " \\\n")
        outputshell.write("-# 100 -n"  + jobIDstring + ".newick\n")
        outputshell.write("\n")
        outputshell.write("mv *.newick " + outputpath + "/\n")
        outputshell.write("rm *_" + str(NumberID) + ".newick.RUN*\n")
    jobcommands.append(jobRunCommand)
    outputfolders.append(outputpath)
with open("CladebreakerCommandList.txt", "w") as OutputListCommands:
    for i in range(len(jobcommands)):
        OutputListCommands.write(jobcommands[i] + "\n")

with open("Cladebreakeroutputfolders.txt", "w") as outputfolderslist:
    for i in range(len(outputfolders)):
        outputfolderslist.write("cp " + outputfolders[i] + "/RAxML_bestTree.CC*  /home/acampbe/DFU/data/WGS_2020/cladebreaker/TreesDFUPatients/\n")
