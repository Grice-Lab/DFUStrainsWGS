import os
import sys
from Bio import SeqIO

#e.g, "/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/IndividualFastas_Hclust"
PathToFastas=sys.argv[1]
# e.g, /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/Hclust_Multifastas/
OutputFolder=sys.argv[2]

ListOfPhageFiles = os.listdir(PathToFastas)

ListOfPhages = set(list(map( lambda x: x.split("_")[0],ListOfPhageFiles)))


for p in ListOfPhages:
    OutputFilePath = os.path.join(OutputFolder, str(p)+"_all.fasta")
    pstring=p+"_"
    SpecificPhageFiles = list(filter(lambda x: pstring in x, ListOfPhageFiles))
    outputobj=open(OutputFilePath, "w")
    for f in SpecificPhageFiles:
        namestring=f.replace(".fasta", "")
        inputfile=os.path.join(PathToFastas, f)
        multifastaobj = list(SeqIO.parse(inputfile, "fasta"))
        for k in range(len(multifastaobj)):
            idstring=namestring + "_"+str(k)
            outputobj.write(">" + idstring)
            outputobj.write('\n')
            outputobj.write(str(multifastaobj[k].seq))
            outputobj.write('\n')
    outputobj.close()
