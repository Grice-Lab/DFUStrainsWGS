# Amy Campbell
# 2023
# Make a multifasta of a representative of every plasmid detected in the dataset
# which will then be concatenated with the phage database
# and then used for blastn-based presence/absence of each gene in the pan-genome

import os
from Bio import SeqIO
import pandas
import sys


# e.g, /home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/PlasmidFastaPaths.csv
PlasmidCSV = sys.argv[1]

AllPlasmidsCSV=sys.argv[2]

# e.g., /home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/PlasmidReps.fasta
outputpathselected = sys.argv[3]

# e.g., "/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/RepMultifastas/"
outputprefixAll = sys.argv[4]

plasmidDF = pandas.read_csv(PlasmidCSV)

AllPlasmidsDF =pandas.read_csv(AllPlasmidsCSV)

for p in set(AllPlasmidsDF['PlasmidID']):
    IndividualPlasmidType = AllPlasmidsDF[AllPlasmidsDF.PlasmidID==p]
    #print(IndividualPlasmidType)
    outputfpath = os.path.join(outputprefixAll, (str(p)+ "_sequences.fasta"))
    #print(outputfpath)
    outputfastaobj = open(outputfpath, "w")
    for f in range(IndividualPlasmidType.shape[0]) :
        inputfastapath = (IndividualPlasmidType.iloc[f])['pathstring']
        inputfastaDORN = (IndividualPlasmidType.iloc[f])['IsolateID']
        fastaobj = list(SeqIO.parse(inputfastapath, "fasta"))
        for k in range(len(fastaobj)):
            idstring= inputfastaDORN+ "_" + str(p) +"_"+str(k)
            seqstring=str(fastaobj[k].seq)
            outputfastaobj.write(">" + idstring)
            outputfastaobj.write('\n')
            outputfastaobj.write(seqstring)
            outputfastaobj.write('\n')

    outputfastaobj.close()

outputfasta = open(outputpath, "w")
for i in range( (plasmidDF.shape[0])):
    plasmidIDobj = ((plasmidDF.iloc[i])["PlasmidID"])
    filepath= ((plasmidDF.iloc[i])["FastaPath"])
    multifastaobj = list(SeqIO.parse(filepath, "fasta"))
    for j in range(len(multifastaobj)):
        newheader = str(plasmidIDobj) + str("_") + str(j)
        sequenceobj = multifastaobj[j]
        outputfasta.write(">" + newheader)
        outputfasta.write('\n')
        outputfasta.write(str(sequenceobj.seq))
        outputfasta.write('\n')



outputfasta.close()

outputbigfasta = open(outputpath, "w")

for i in range( (plasmidDF.shape[0])):
    plasmidIDobj = ((plasmidDF.iloc[i])["PlasmidID"])
    filepath= ((plasmidDF.iloc[i])["FastaPath"])
    multifastaobj = list(SeqIO.parse(filepath, "fasta"))
    for j in range(len(multifastaobj)):
        newheader = str(plasmidIDobj) + str("_") + str(j)
        sequenceobj = multifastaobj[j]
        outputfasta.write(">" + newheader)
        outputfasta.write('\n')
        outputfasta.write(str(sequenceobj.seq))
        outputfasta.write('\n')



outputfasta.close()
