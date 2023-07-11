# Amy Campbell
# 2023
# Parsing and renaming output of CD-hit-est clustering 
# of all complete phages identified from the 220 S. aureus genomes

# Takes in (positionally):
    # 1. Path to .clstr file from CD-hit-est
    #    (e.g., /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/ClusteredPhages90.clustr)
    # 2. Path to .fasta file from CD-hit-est
    #    (e.g., /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/ClusteredPhages90)
    # 3. Path to folder containing full contigs files for all genomes (suffix "_Final.fasta")
    #   (e.g, home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates)
    # 4. CSV file and path to output presence/absence of each phage in each genome based only on CD-hit-est
    #    (e.g., /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/PhagePresenceAbsence.csv)
    # 5. Fasta file to write the sequences with cleaned up IDs
    #    /home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/ClusteredPhagesRenamed.fasta)

# Outputs:
    # 1. CSV where rows are genome IDs (e.g. DORN76) and columns are phages (Phage0..Phage68), values aare presence(1) or absence (0)
    #   based only on whether their own sequence clustered into the CD-hit-est cluster
    # 2. Multifasta where ids are Phage0 .... Phage68 instead of their current ids

import os
from Bio import SeqIO
import pandas
import sys

ClusteredPhages = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/ClusteredPhages90.clstr"
PhagesFasta = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/ClusteredPhages90"
FinalSetIsolates="/Users/amycampbell/Documents/FinalSetDFUIsolates"
OutputCSV = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/PhagePresenceAbsence.csv"
OutputFasta = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/ClusteredPhagesRenamed.fasta"

allIsolates=os.listdir(FinalSetIsolates)

ClusteredPhageString = open(ClusteredPhages, "r").read()
ClusteredPhageList = ClusteredPhageString.split("\n>")


Multifasta = list(SeqIO.parse(PhagesFasta, "fasta"))



phagelist = (list(map(lambda x: "Phage" + str(x), (list(range(69))))))
DORNlist = list(map(lambda x: str.replace(x, "_Final.fasta", ""), allIsolates))
DORNlist = list(filter(lambda x: "DORN" in x, DORNlist)) # remove the .DSstore file (classic)


# Make an empty dataframe where columns are new phage names, rows are isolate names
# Will be updated to 1 if an isolate contains a given phage cluster
# Also make a dictionary where keys will be new phage IDs (Phage1...Phage68) and values
# will be the sequences associated with those phages
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
            genome = ((((rowitem.split(">"))[1].split("|"))[0]).split("_"))[0]
            PresenceAbsenceDF.at[genome,phageitemstring] = 1

PresenceAbsenceDF.to_csv(OutputCSV)


with open(OutputFasta, "w") as outputfasta:
    for k in NewMultifastaDict.keys():
        outputfasta.write(">" + str(k))
        outputfasta.write('\n')
        outputfasta.write(str(NewMultifastaDict[k]))
        outputfasta.write('\n')
