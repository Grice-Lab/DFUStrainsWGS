# Amy Campbell
# 2022
# Looking for NS mutations in rsbU AA sequences in the DFU collection
# IDC about efficiency here, just correctness

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import os


def IDmismatches(characterlistQ, characterlistR, Qname):
    # characterlistR character list like ["K", "L", "G", "M"] of reference rsbU sequence (USA300/LAC)
    # characterlistQ character list of query rsbU sequence (whatever mismatching version)
    # Qname is the genome name (DORN)
    # Returns list of lists [[]] where each inner list is an instance of non-matching:
    # [Qname, NonmatchPosition, ReferenceLetter, QueryLetter] where position is 1-indexed (sry python)
    result=[]
    for i in range(len(characterlistQ)):
        if characterlistQ[i] != characterlistR[i]:
            result.append([Qname, i+1, characterlistR[i], characterlistQ[i]])
    return result

referenceRsbU = "/Users/amycampbell/Documents/GriceLabGit/DFUStrainsWGS/Association/rsbU_all/FPR3757rsbU.faa"
AllFAAs = "/Users/amycampbell/Documents/DFUData/AllFAAs"

recordRsbURef = list(SeqIO.parse(referenceRsbU, "fasta"))[0]

referenceSeq = recordRsbURef.seq
faaPrefix="/Users/amycampbell/Documents/DFUData/AllFAAs/"

# Grab rsbU for each genome
###########################

FileList = (os.listdir(AllFAAs))

outputseqs = []
for f in FileList:
    justName = f.replace(".faa", "")

    records = SeqIO.to_dict(SeqIO.parse(faaPrefix+f, "fasta"), key_function=lambda rec : rec.description)
    records_containing_rsbU =  [desc for desc in records.keys() if "PP2C family protein-serine/threonine phosphatase" in desc]

    if len(records_containing_rsbU) > 1:
        print("Oops more than 1 rsbU eek")
    else:
        if len(records_containing_rsbU) < 1:
            print(justName)
            print("has none!")
        else:
            sequence_rsbu = records[records_containing_rsbU[0]].seq
            #print(len(sequence_rsbu))
            id_rsbu = justName + "_rsbU"
            newseqrec = SeqIO.SeqRecord(id=id_rsbu,seq=sequence_rsbu )
            outputseqs.append(newseqrec)

# Write out .faa file with all the rsbU sequences
#################################################

with open("AllRsbUs.faa", "w") as outputfasta:
    SeqIO.write(outputseqs,outputfasta, "fasta")

# For each record of rsbU we have,

##########################################
Identicals=[]
NonIdenticals=[]

for record in outputseqs:

    if(str(record.seq)==str(referenceSeq)):
        Identicals.append(record)
    else:
        NonIdenticals.append(record)

characterlistReference = list(referenceSeq)
print(characterlistReference)
BigStupidList = []
for record in NonIdenticals:
    characterlistQuery = list(record.seq)
    QName=record.id
    QName = QName.replace("_rsbU", "")
    MismatchList = IDmismatches(characterlistQuery, characterlistReference,QName)
    BigStupidList.extend(MismatchList)


with open("RsbUMismatches.csv", "w") as outputcsv:
    outputcsv.write("Genome, NonmatchPos, ReferenceAA, QueryAA\n")
    for l in range(len(BigStupidList)):
        outputstring= ','.join(map(str, BigStupidList[l]))
        outputcsv.write(outputstring)
        outputcsv.write("\n")
