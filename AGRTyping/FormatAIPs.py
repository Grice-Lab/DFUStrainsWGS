# Amy Campbell
# 02/2020

# Takes in the AGRs_CloseScores.csv,
# which lists genomes which contain both AIP sequences
# as well as their protein IDs (from prokka output)

# For each genome in 'Genome' column:
    # Opens the corresponding GENOMENAME.tbl file from
    # /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/GENOMENAME/___
    # Splits it up by occurences of ">Feature" (demarcating contigs)

    # Set Assignment1, Assignment2, Protein1, Protein2, GenomeName
    # Set 'foundcounts' to 0
    # Loop 1: getting protein IDs
    ###################
    # While foundcounts <2:
    #   go through items in this list.
    #   Check for substring BestProteinID.
    #       If found, make a dictionary entry where the key is GENOMENAME__BestProteinID__Assignment,
    #       item is an array of the protein ids of everything in that feature (delimited by 'locus_tag       ')
    #       Increment foundcounts by 1.
    #   Check for substring NextBestProteinID
    #       If found,  make a dictionary entry where the key is GENOMENAME__NextBestProteinID__NextBestScoringAssignment,
    #       item is an array of the protein ids of everything in that features (delimited by 'locus_tag       ')
    #       Increment foundcounts by 1.
    #
    # Close the GENOMENAME.tbl file.
    ###################


#######################
import pandas
import os
from Bio import SeqIO

def ProteinIDList(FeatureString):
    #update
    FeatureStringList = FeatureString.split('\tlocus_tag\t')
    #print(FeatureStringList)
    returnList = list(map(lambda substr: (substr.split())[0],FeatureStringList[1:]))
    #print(returnList)
    return returnList

AIP_IDs_FilePath ="/home/acampbe/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/AGRs_CloseScores.csv"
PathPrefix="/home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files"

AIPTable = pandas.read_csv(AIP_IDs_FilePath, header=0)
ProteinDictionary=dict()
genomes= list(AIPTable['Genome'])
PRoteinIDDict=dict()
#PRoteinIDDict
PRoteinAssignmentDict = dict()
for index, GenomeRow in AIPTable.iterrows():
    GenomeName = str(GenomeRow['Genome'])
    Assignment1 = str(GenomeRow['Assignment'])
    Assignment2 = str(GenomeRow['NextBestScoringAssignment'])
    Protein1 = str(GenomeRow['BestProteinID'])
    Protein2 = str(GenomeRow['NextBestProteinID'])
    print(GenomeName)
    PRoteinIDDict[GenomeName] = [Protein1, Protein2]
    PRoteinAssignmentDict[GenomeName] = [Assignment1, Assignment2]
    tblpath=os.path.join(PathPrefix, GenomeName)
    tblfile=os.path.join(tblpath, str(GenomeName + ".tbl"))
    tblObject = open(tblfile, "r")
    tblString = tblObject.read()
    tblObject.close()

    tblStringList = tblString.split(">Feature")

    foundcounts=0
    featindex=0
    while foundcounts < 2:
        string_item = str(tblStringList[featindex])
        if Protein1 in string_item:
            ProteinDictionary[str(GenomeName + "__" + Protein1 + "__" + Assignment1)] = ProteinIDList(string_item)
            foundcounts = foundcounts + 1
        elif Protein2 in string_item:
            ProteinDictionary[str(GenomeName + "__" + Protein2 + "__" + Assignment2)] = ProteinIDList(string_item)
            foundcounts = foundcounts + 1
        featindex = featindex + 1


for genome in genomes:
    genome = str(genome)
    folderpath = os.path.join(PathPrefix, genome)
    faapath = os.path.join(folderpath, str(genome + ".faa"))
    outputpath= os.path.join(folderpath, str(genome + "_AGRContigs.faa" ))
    fastaObject= SeqIO.parse(faapath, "fasta")
    outputfile = open(outputpath, "w")
    for record in fastaObject:
        infolist=str(record.id).split(" ")
        if infolist[0] in ProteinDictionary[str(genome + "__" + (PRoteinIDDict[genome])[0] + "__" +  (PRoteinAssignmentDict[genome])[0])]:
            record.id = str(genome + "_" + (PRoteinAssignmentDict[genome])[0] + "_" + record.id)
            SeqIO.write(record, outputfile, 'fasta')
        elif infolist[0] in ProteinDictionary[str(genome + "__" + (PRoteinIDDict[genome])[1] + "__" +  (PRoteinAssignmentDict[genome])[1])]:
            record.id = str(genome + "_" + (PRoteinAssignmentDict[genome])[1] + "_" + record.id)
            SeqIO.write(record, outputfile, 'fasta')
    outputfile.close()




#
#fastaObject= SeqIO.parse(gffpath, "fasta")
#unfound=len(combinedlist)
#while unfound > 0:







    #
