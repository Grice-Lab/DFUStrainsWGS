# Amy Campbell
# Make a mapping of the roary pan genome results with vs. without 2149 roary so you can plot the 2149-present pan-genome but use the
# Blast2GO classifications, which were a nightmare to obtain and were done using the non-2149 version

from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

#from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas
import subprocess
import os
import inspect
import sys

CorrectPanGenome = "/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/pan_genome_reference_DORN2149.fa"
OldPanGenome = "/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022_old/pan_genome_reference_roary2022_old.fa"
OutputBlastDB = "/Users/amycampbell/Documents/DataInputGithub/data/AccessoryGenomeMap.xml"


NewPanGenomeSeqs = SeqIO.to_dict(SeqIO.parse(open(CorrectPanGenome),'fasta'))


listRecords = []

# Make a dictionary of roary orthol

# Database is the old pan-genome
################################
#commandMakeDB = "makeblastdb -in " + OldPanGenome + " -dbtype nucl -input_type fasta -out " + OutputBlastDB
#subprocess.run(commandMakeDB, shell=True)

# Query is new(2149-inclusive) pan-genome
#########################################
#commandRunBlast ="blastn -query " + CorrectPanGenome + " -db " + OutputBlastDB + " -outfmt 5 -out Old_New_Pangenome_blasthits.xml"
#subprocess.run(commandRunBlast, shell=True)

result_handle = open("Old_New_Pangenome_blasthits.xml", 'r')
blast_records = NCBIXML.parse(result_handle)
outputfasta=open("/Users/amycampbell/Documents/DataInputGithub/data/NeedsBlast2GOannotations.fa", "w")
outputCSV = open("/Users/amycampbell/Documents/DataInputGithub/data/Map_New_To_Old_Pangenome.csv", "w")
# Keys are 'new', values are 'old'

MappingDict = dict()

for blast_record in blast_records:
    if len(blast_record.alignments) > 0 :
        TopIdentity = (blast_record.alignments[0].hsps[0].identities)/blast_record.query_letters
        #print(TopIdentity)
        #print(blast_record.query.split(" "))
        OldName = ((blast_record.alignments[0].hit_def).split(" "))[0]
        NewName = (blast_record.query.split(" "))[0]
        MappingDict[NewName] = OldName


        if(TopIdentity < .7):
            NewName = (blast_record.query.split(" "))[0]
            print(NewName)
            print("best hit has " + str(round(100*TopIdentity, 3)) + "% identity" )
            MappingDict[NewName] = "New"

    else:
        NewName = (blast_record.query.split(" "))[0]
        SeqIO.write(NewPanGenomeSeqs[NewName], outputfasta, "fasta")
        MappingDict[NewName] = "New"

outputCSV.write("NewpangenomeID,OldPangenomeID\n")
for k in MappingDict.keys():
    outputCSV.write(k + ","+ MappingDict[k] + "\n")
outputCSV.close()
outputfasta.close()
