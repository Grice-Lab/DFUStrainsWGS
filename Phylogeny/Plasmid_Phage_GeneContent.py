from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas
import subprocess
import os
import inspect
import sys
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline


# pangenomeReferencePath = sys.argv[1]
# PhagePlasmidFasta = sys.argv[2]
# blastDBpath = sys.argv[3]
# outputDFpath = sys.argv[4]


pangenomeReferencePath="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/pan_genome_reference_DORN2149.fa"
PhagePlasmidFasta = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/PhagesPlasmids.fasta"
blastoutputpath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/plasmid_phage_blast/"
outputDFpath = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/GenePresence_ByPhage_Plasmid.csv"
def getgeneName(seqobj):
    return( ((seqobj.description).split(" "))[1])

def getHGTnames(seqobj):
    return(((seqobj.id).split("_"))[0])

def getgeneLength(seqobj):
    return(len(str(seqobj.seq)))
# iterate through the pan_genome

# 1. Make blastdb
###################

blastDBpath = os.path.join(blastoutputpath, "phagePlasmidDB")
commandMakeBlastDB = "makeblastdb -in " + PhagePlasmidFasta + " -dbtype nucl -input_type fasta -out " + blastDBpath
#subprocess.run(commandMakeBlastDB, shell=True)


# 2. Read in the multifasta of pan-genome .fna sequences
#########################################################
pangenomeseqs = list(SeqIO.parse(open(pangenomeReferencePath), 'fasta'))
PhagePlasmidSeqs = list(SeqIO.parse(open(PhagePlasmidFasta), 'fasta'))

genenames = list(map(getgeneName, pangenomeseqs))
genelengths = list(map(getgeneLength, pangenomeseqs))
plasmidPhagenames = set(list(map(getHGTnames, PhagePlasmidSeqs)))
PresenceAbsenceDF = pandas.DataFrame(0, index=genenames, columns=plasmidPhagenames)

# Just so I'm not carrying this gigantic list around
PhagePlasmidSeqs = list()

with open(pangenomeReferencePath.replace(".fa", "_simplified.fasta"), "w" ) as simplifiedfasta:
    for item in pangenomeseqs:
        gene = str(((item.description).split(" "))[1])
        simplifiedfasta.write(">"+gene)
        simplifiedfasta.write("\n")
        simplifiedfasta.write(str(item.seq))
        simplifiedfasta.write("\n")

# Run blastn
############
inputQueryseqs = pangenomeReferencePath.replace(".fa", "_simplified.fasta")
blastn_output = os.path.join(blastoutputpath, "pangenome_blast_HGT.tab")
commandblastn = "blastn -query " + inputQueryseqs + " -db "+ blastDBpath +" -outfmt 6 -out " + blastn_output
#subprocess.run(commandblastn, shell=True)

blastoutput_df = pandas.read_csv(blastn_output, sep="\t")

blastoutput_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

blastoutput_df['HGTid']=  blastoutput_df['sseqid'].apply(lambda x: x.split("_")[0])# (map(blastoutput_df.sseqid, lambda x: x.split("_")[0]))


for i in range(len(genenames)):
    glength = genelengths[i]
    gname = genenames[i]
    subsetdf = blastoutput_df[blastoutput_df['qseqid']==gname]
    subsetdf['pctcovered'] = subsetdf['length'].apply(lambda x: 100*(x/glength)) #df.apply(lambda x: (x['sseqid'].split("_"))[0])
    subsetdf = subsetdf[ ((subsetdf['pident']> 90) & (subsetdf['pctcovered']> 90) ) ]
    if (subsetdf.shape[0]) > 0:
        PresenceAbsenceDF.loc[gname,subsetdf['HGTid'] ] = 1

PresenceAbsenceDF.to_csv(outputDFpath)

    #        PresenceAbsenceDF.loc[gname, ]
# # read blastn output
# result_handle = open(blastn_output, 'r')
# blast_records = NCBIXML.parse(result_handle)
# i=1
# for recordblast in blast_records:
#     geneitem = (recordblast.query.split(" "))[0]
#     if geneitem=="sak":
#         lengthquery = recordblast.query_letters
#         for alignment in recordblast.alignments:
#             print(alignment.hit_def)
#             print(alignment.hsps)
