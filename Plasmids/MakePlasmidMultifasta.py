# Amy Campbell
# 2023
# Make a multifasta of a representative of every plasmid detected in the dataset
# which will then be concatenated with the phage database
# and then used for blastn-based presence/absence of each gene in the pan-genome

import os
from Bio import SeqIO
import pandas
import sys

PlasmidCSV = sys.argv[1]
outputpath = sys.argv[2]


plasmidDF = pandas.read_csv(PlasmidCSV)


outputfasta = open(outputpath, "w")
for i in range( (plasmidDF.shape[0])):
    plasmidIDobj = ((plasmidDF.iloc[i])["PlasmidID"])
    filepath= ((plasmidDF.iloc[i])["FastaPath"])
    multifastaobj = list(SeqIO.parse(filepath, "fasta"))
    for(j in range(len(multifastaobj))):
        newheader = str(plasmidIDobj) + str("_") + str(j)
        sequenceobj = multifastaobj[j]
        outputfasta.write(">" + newheader)
        outputfasta.write('\n')
        outputfasta.write(str(sequenceobj.seq))
        outputfasta.write('\n')



outputfasta.close()
