# Amy Campbell
# July 2023
# Make multifastas of Phages
# Takes in:
#   folder path containing DORN933, DORN1088,...etc folders of checkV output files
#   folder path for output


import os
from Bio import SeqIO
import pandas
import sys


InputPath= sys.argv[1]
OutputPath=sys.argv[2]

folderlist = os.listdir(InputPath)

MultifastaDict = dict()
PDlist = []

for foldername in folderlist:
    FullPath = os.path.join(InputPath,foldername )
    Genome=foldername
    QualitySummary = pandas.read_table(os.path.join(FullPath, "quality_summary.tsv"))
    #print(QualitySummary['checkv_quality'])
    QualitySummary = QualitySummary[(QualitySummary['completeness']==100.0) & ((QualitySummary['checkv_quality']=="High-quality") | (QualitySummary['checkv_quality']=="Complete"))]

    if(QualitySummary.empty):
        print("No complete phages in "+Genome )
    else:
        viraldict=SeqIO.to_dict(SeqIO.parse(os.path.join(FullPath, "viruses.fna"), "fasta"))
        proviraldict=SeqIO.to_dict(SeqIO.parse(os.path.join(FullPath, "proviruses.fna"), "fasta"))
        for index in range( QualitySummary.shape[0] ):
            contigid=QualitySummary.iloc[index].contig_id
            provirus_status=QualitySummary.iloc[index].provirus
            if provirus_status=="Yes":
                contigid=str(contigid) + "_1"
                sequence = proviraldict[contigid].seq
                newcontigid=str(Genome) + "_" + contigid
                MultifastaDict[newcontigid] = sequence
                lengthobj=QualitySummary.iloc[index].proviral_length
                Type="provirus"

            else:
                sequence = viraldict[contigid].seq
                newcontigid=str(Genome) + "_" + contigid
                MultifastaDict[newcontigid] = sequence
                lengthobj=QualitySummary.iloc[index].contig_length
                Type="virus"

            PDlist.append([Genome, newcontigid, Type, lengthobj])

PD_summary = pandas.DataFrame(PDlist, columns=["Genome", "SequenceName", "SequenceType", "Length"])
PD_summary.to_csv(os.path.join(OutputPath, "SummaryViralSequences.csv"))

with open(os.path.join(OutputPath, "AllPhages.fasta"), "w") as outputfasta:
    for k in MultifastaDict.keys():
        outputfasta.write(">" + str(k))
        outputfasta.write('\n')
        outputfasta.write(str(MultifastaDict[k]))
        outputfasta.write('\n')

















#
