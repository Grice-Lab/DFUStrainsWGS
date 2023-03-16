# Amy Campbell
# March 2023
# Testing the 'top hits' of WhatsGNU output  for uniqueness within the genome to reduce ambiguous mapping issues

# Takes in:
# - some genome's WhatsGNU WhatsGNU report (which includes mappings of annotation protein IDs --> WG ortholog groups)
# - List of the whatsGNU ortholog groups to test
# - Roary 'clustered proteins' file (which shows which Roary ortholog IDs include which protein IDs)
# - .faa of all the annotated proteins in the genome to act as the blast database

# Then:
# Makes a dictionary mapping the WG orthologs to test (keys) to the genome's protein IDs (values)
# reads in the .faa file of the genome as a
# Iterates through that dictionary's keys, and for each one:
#   Gets the protein ID corresponding to it
#   Does a blast search of that protein ID's sequence against all the protein sequences for the genome
#   If there's no good hits to this OTHER than itself, INCLUDE the protein by:
#       Adding the protein ID, its corresponding WG ortholog group ID,and its roary ortholog group as a row to a data frame
#   If there are good hits, exclude


# Outputs:

from Bio.Blast.Applications import NcbiblastpCommandline
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

pangenomepath= "/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/pan_genome_reference.fa"
orthologlistpath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/Patient176_ortholog_list.txt"
clusteredproteinspath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/clustered_proteins"
WhatsGreport="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/DORN1646_WhatsGNU_report.txt"
GenomeFAApath = "/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/DORN1646.faa"
GenomeDB = os.path.dirname(orthologlistpath)
Genomestring=os.path.basename(GenomeFAApath).replace(".faa", "")
GenomeDBpath = os.path.join(GenomeDB, Genomestring + "_blast" )
print(GenomeDB)
#
# # Make a blast database of the genome's faa
commandMakeDB = "makeblastdb -in " + GenomeFAApath + " -dbtype prot -input_type fasta -out " + GenomeDBpath
subprocess.run(commandMakeDB, shell=True)

# get a list of ortholog IDs from whatsgnu to include
# (output by WhatsGNUoutput.R)
#####################################################
readorthologs  = open(orthologlistpath,"r" ).read().split('\n')
readorthologs.remove("")
# Load the WhatsGNU report for a genome
WGtable = pandas.read_csv(WhatsGreport,sep='\t')

# Just get the protein id out of there
######################################
WGtable["protID"] = WGtable['protein'].str.split(" ").str[0]

# keys are the WG ortholog groups, values are the annotated protein IDs
########################################################################
MappingWG_protID= dict(zip(WGtable.ortholog_group, WGtable.protID))

# filter to just those which are in the list of 'differentiating' orthologs
###########################################################################
MappingWG_protID_filtered = dict((k, MappingWG_protID[k]) for k in readorthologs)

# Make a dictionary of all the proteins and their AA sequences we're considering using
######################################################################################
GenomeSeqHitsDB = dict()
GenomesProteins = SeqIO.parse(open(GenomeFAApath),'fasta')
listRecords = []

# Make a dictionary of roary ortholog names to the proteins they contain
clusterdict = dict()
clusteredproteins = open(clusteredproteinspath,"r").read().split('\n')
clusteredproteins.remove('')

for c in clusteredproteins:
    littlelist=c.split(": ")
    orthname = littlelist[0]
    proteins=littlelist[1].split('\t')
    clusterdict[orthname] = proteins

Hitsfasta = Genomestring + "_hits.faa"
for p in GenomesProteins:
    if p.id in MappingWG_protID_filtered.values():
        GenomeSeqHitsDB[p.id] = p.seq
        listRecords.append(p)

# Write it out as <GenomeName>_hits.faa
SeqIO.write(listRecords, Hitsfasta, "fasta")

# Run blastP for each query sequence in there against the database of all the proteins in the genome
####################################################################################################
commandRunBlastP ="blastp -query " + Hitsfasta + " -db " + GenomeDBpath + " -outfmt 5 -out "+ Genomestring +  "_blasthits.xml"
print(commandRunBlastP)
subprocess.run(commandRunBlastP, shell=True)


result_handle = open(Genomestring +  "_blasthits.xml", 'r')
blast_records = NCBIXML.parse(result_handle)

ProteinsKeep = []
for blast_record in blast_records:
    protein_query_name = (blast_record.query.split(" "))[0]

    alignmentnum=0
    # % identity of first hit for the query (which should just be itself)
    FirstIdentity = (blast_record.alignments[0].hsps[0].identities)/blast_record.query_letters
    if len(blast_record.alignments)>1:
        # % identity of second best hit for the query
        SecondIdentity = (blast_record.alignments[1].hsps[0].identities)/blast_record.query_letters

        # Don't want to use genes that have >70% Protein identity to anything other than itself
        if SecondIdentity > .6:
            Keep=False
        else:
            Keep=True
    # if the only alignment a gene has is to itself then include it
    else:
        Keep=True

    # only keep a protein if it passed the uniqueness test
    if Keep==True:
        ProteinsKeep.append(protein_query_name)

outputfile= open(os.path.join(GenomeDB,os.path.basename(GenomeDB) + "_unique.txt"), "w")

# Look up which roary ortholog this corresponds to
clusterkeys = list(clusterdict.keys())
print("# Proteins kept: " +str(len(ProteinsKeep)) )
for p in ProteinsKeep:
    found=False
    it=0
    while(found==False):
        #print(clusterdict[clusterkeys[it]])
        if (p in list(clusterdict[clusterkeys[it]])):
            RoaryOrth = (clusterkeys[it])
            outputfile.write(RoaryOrth+'\n')
            found=True

        if (it>len(clusterdict)) & (found==False):
            print("error")
            found=True

        it=it+1
