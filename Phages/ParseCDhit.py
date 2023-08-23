# Amy Campbell
# 2023
# Parsing and renaming output of CD-hit-est clustering
# of all complete phages identified from the 220 S. aureus genomes

# Takes in (positionally):
    # 1. Path to .clstr file from CD-hit-est
    #    (e.g., /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/ClusteredPhages80.clustr)
    # 2. Path to .fasta file from CD-hit-est
    #    (e.g., /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/ClusteredPhages80)
    # 3. Path to folder containing full contigs files for all genomes (suffix "_Final.fasta")
    #   (e.g, home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates)
    # 4. CSV file and path to output presence/absence of each phage in each genome based only on CD-hit-est
    #    (e.g., /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/PhagePresenceAbsence.csv)
    # 5. Fasta file to write the sequences with cleaned up IDs
    #    /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/ClusteredPhagesRenamed.fasta)
    # 6. Folder to write individual fastas to (named like)
    #  e.g., /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/IndividualFastas/
    # 7. Folder containing all the phages input into the clustering (/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/AllPhages.fasta)
# Outputs:
    # 1. CSV where rows are genome IDs (e.g. DORN76) and columns are phages (Phage0..Phage68), values aare presence(1) or absence (0)
    #   based only on whether their own sequence clustered into the CD-hit-est cluster
    # 2. Multifasta where ids are Phage0 .... Phage68 instead of their current ids

import os
from Bio import SeqIO
import pandas
import sys

ClusteredPhages = sys.argv[1]
PhagesFasta = sys.argv[2]
FinalSetIsolates=sys.argv[3]
OutputCSV = sys.argv[4]
OutputFasta = sys.argv[5]
OutputFolderFastas=sys.argv[6]
InputFullPhageList=sys.argv[7]


allIsolates=os.listdir(FinalSetIsolates)

ClusteredPhageString = open(ClusteredPhages, "r").read()
ClusteredPhageList = ClusteredPhageString.split("\n>")

# List of all representative sequences
Multifasta = list(SeqIO.parse(PhagesFasta, "fasta"))

phagelist = (list(map(lambda x: "Phage" + str(x), (list(range(len(Multifasta)))))))
DORNlist = list(map(lambda x: str.replace(x, "_Final.fasta", ""), allIsolates))
DORNlist = list(filter(lambda x: "DORN" in x, DORNlist)) # remove the .DSstore file (classic)

allphagesdict = dict()
AllUnclusteredPhages =  list(SeqIO.parse(InputFullPhageList, "fasta"))
for j in range(len(AllUnclusteredPhages)):
    fastitem = AllUnclusteredPhages[j]
    itemID=(AllUnclusteredPhages.id).split("|")[0]
    allphagesdict[itemID] = fastitem.seq
AllUnclusteredPhages=""
# Make an empty dataframe where columns are new phage names, rows are isolate names
# Will be updated to 1 if an isolate contains a given phage cluster
# Also make a dictionary where keys will be new phage IDs (Phage1...Phage52) and values
# will be the representative sequences associated with those phages
###################################################################################
PresenceAbsenceDF = pandas.DataFrame(0, index=DORNlist, columns=phagelist)
NewMultifastaDict = dict()

for s in range(len(Multifasta)):
    FastaItem = Multifasta[s]
    IdentifyingString = (FastaItem.id).split(":")[0]
    IdentifyingString = (IdentifyingString + ":")
    MatchedItems = [item for item in ClusteredPhageList if IdentifyingString in item]
    ClusterItem = MatchedItems[0]
    ClusterItemList = ClusterItem.split('\n')
    ClusterNum = ClusterItemList[0].split(" ")[1]

    phageitemstring = "Phage"+ str(ClusterNum)

    # Make new entry for this phage in the new multifasta dict
    NewMultifastaDict[phageitemstring] = FastaItem.seq

    # Mark this phage as being present in the genomes it's present in
    for i in range(1,len(ClusterItemList)):
        rowitem = ClusterItemList[i]
        if rowitem != "":
            phageinstanceSearchString =  str(((rowitem.split(">"))[1].split("|"))[0])

            genome = ((((rowitem.split(">"))[1].split("|"))[0]).split("_"))[0]

            PresenceAbsenceDF.at[genome,phageitemstring] = 1

            sequenceitem  = allphagesdict[phageinstanceSearchString]
            newid = phageitemstring+"_" + phageinstanceSearchString
            outputlittlefasta=open(OutputFolderFastas + newid + ".fasta", "w")
            outputlittlefasta.write(">"  + newid)
            outputlittlefasta.write('\n')
            outputlittlefasta.write(str(sequenceitem))
            outputlittlefasta.close()
PresenceAbsenceDF.to_csv(OutputCSV)


with open(OutputFasta, "w") as outputfasta:
    for k in NewMultifastaDict.keys():
        outputfasta.write(">" + str(k))
        outputfasta.write('\n')
        outputfasta.write(str(NewMultifastaDict[k]))
        outputfasta.write('\n')
