# Amy Campbell
# Extract yabJ-spoVG genes from the
# gene-level alignment from Roary

import os
from Bio import SeqIO
import pandas
import sys

spoVGpath = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/SpoVGmappings.csv"
yabJpath = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/YabJmappings.csv"

yabJpangenomePath="/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/extdb_pgaptmp_001576.fa.aln"
spoVG_pangenomePath = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/spoVG.fa.aln"

yabJmap = pandas.read_csv(yabJpath)
spoVGmap =  pandas.read_csv(spoVGpath)

yabJoutputfasta="/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/YabJSequences.fasta"
spoVGoutputfasta="/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/SpoVGSequences.fasta"

yabJsequences = SeqIO.to_dict(SeqIO.parse(yabJpangenomePath, "fasta"))
spoVGsequences = SeqIO.to_dict(SeqIO.parse(spoVG_pangenomePath, "fasta"))


yabJoutput = open(yabJoutputfasta, "w")
for rownum in range(yabJmap.shape[0]):
    genomestring = (yabJmap.iloc[rownum].Genome)
    keystring=(yabJmap.iloc[rownum].Mapping)
    seqstring = str(yabJsequences[keystring].seq)
    seqstring = seqstring.replace("-","")
    yabJoutput.write(">" + genomestring +"\n")
    yabJoutput.write(seqstring + "\n")

yabJoutput.close()


spoVGoutput = open(spoVGoutputfasta, "w")
for rownum in range(spoVGmap.shape[0]):
    genomestring = (spoVGmap.iloc[rownum].Genome)
    keystring=(spoVGmap.iloc[rownum].Mapping)
    seqstring = str(spoVGsequences[keystring].seq)
    seqstring = seqstring.replace("-","")
    spoVGoutput.write(">" + genomestring +"\n")
    spoVGoutput.write(seqstring + "\n")

spoVGoutput.close()
