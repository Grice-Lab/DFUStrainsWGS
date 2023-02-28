# Amy Campbell
# Feb 2023
# Make a shellscript for each patient with a CSV in input csvs folder for cladebreaker
import os

PathToInputs="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/InputCSVs/"
configstring="/home/acampbe/DFU/data/WGS_2020/cladebreaker/penn_cluster.config"


listInputs = os.listdir(PathToInputs)
print(listInputs)

for l in range(len(listInputs)):
    inputcsvpath = PathToInputs + str(listInputs[l])
    jobIDstring=listInputs[l].replace("_input.csv", "")

print(jobIDstring)
