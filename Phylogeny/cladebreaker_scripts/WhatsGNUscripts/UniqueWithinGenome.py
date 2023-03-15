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
# pangenomepath= "/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/pan_genome_reference.fa"
# orthologlistpath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/Patient176_ortholog_list.txt"
# clusteredproteinspath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/clustered_proteins"
# WhatsGreport="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/DORN1646_WhatsGNU_report.txt"
# GenomeFAApath = "/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/DORN1646.faa"
# GenomeDB = os.path.dirname(orthologlistpath)
# Genomestring=os.path.basename(GenomeFAApath).replace(".faa", "")
# GenomeDBpath = os.path.join(GenomeDB, Genomestring + "_blast" )
#
#
# # Make a blast database of the genome's faa
# commandMakeDB = "makeblastdb -in " + GenomeFAApath + " -dbtype prot -input_type fasta -out " + GenomeDBpath
# subprocess.run(commandMakeDB, shell=True)
#
# # get a list of ortholog IDs from whatsgnu to include
# # (output by WhatsGNUoutput.R)
# #####################################################
# readorthologs  = open(orthologlistpath,"r" ).read().split('\n')
# readorthologs.remove("")
# # Load the WhatsGNU report for a genome
# WGtable = pandas.read_csv(WhatsGreport,sep='\t')
#
# # Just get the protein id out of there
# ######################################
# WGtable["protID"] = WGtable['protein'].str.split(" ").str[0]
#
# # keys are the WG ortholog groups, values are the annotated protein IDs
# ########################################################################
# MappingWG_protID= dict(zip(WGtable.ortholog_group, WGtable.protID))
# # filter to just those which are in the list of 'differentiating' orthologs
# ###########################################################################
# MappingWG_protID_filtered = dict((k, MappingWG_protID[k]) for k in readorthologs)
#
# # Make a dictionary of all the proteins and their AA sequences we're considering using
# ######################################################################################
# GenomeSeqHitsDB = dict()
# GenomesProteins = SeqIO.parse(open(GenomeFAApath),'fasta')
# listRecords = []
#
# Hitsfasta = Genomestring + "_hits.faa"
# for p in GenomesProteins:
#     if p.id in MappingWG_protID_filtered.values():
#         GenomeSeqHitsDB[p.id] = p.seq
#         listRecords.append(p)
#
# # Write it out as <GenomeName>_hits.faa
# SeqIO.write(listRecords, Hitsfasta, "fasta")
#
# # Run blastP for each query sequence in there against the database of all the proteins in the genome
# ####################################################################################################
# commandRunBlastP ="blastp -query " + Hitsfasta + " -db " + GenomeDBpath + " -outfmt 5 -out "+ Genomestring +  "_blasthits.xml"
# print(commandRunBlastP)
# subprocess.run(commandRunBlastP, shell=True)
#
# clusteredproteins = open(clusteredproteinspath,"r").read().split('\n')
#
#
#
# proteinstring="KEKMFCID_00297"
# GenomeSeqHitsDB[proteinstring]

result_handle = open("DORN1646_blasthits.xml", 'r')
blast_records = NCBIXML.parse(result_handle)


# look to see how many hits there are that have >.5*bit score of the
#
recordnum=0
for blast_record in blast_records:
    protein_query_name = (blast_record.query.split(" "))[0]

    alignmentnum=0
    FirstIdentity = (blast_record.alignments[0].hsps[0].identities)/blast_record.query_letters
    if len(blast_record.alignments)>1:
        SecondIdentity = (blast_record.alignments[1].hsps[0].identities)/blast_record.query_letters

        if SecondIdentity > .8:
            print(blast_record.query)
            print(blast_record.alignments[1].hsps[0].sbjct)
    else:
        print("Just one alignment")
    # if(recordnum==0):
    #     #print(blast_record.query_letters)
    #     print(FirstIdentity)
    #     print(SecondIdentity)
    #
    # for align in blast_record.alignments:
    #     if alignmentnum==0 and recordnum==0:
    #         print("")#print("identity " + str(align.hsps[0].identities))
    #
    #
    # recordnum=recordnum+1
        # for hsp in align.hsps:
        #     if
        #
        #     alignmentnum=alignmentnum+1

    # Do something with blast_record
#GenomesProteins
###################################
#for orthgroups in MappingWG_protID.keys():

# alignmentnum=alignmentnnum+1
# for hsp in align.hsps:
#     if alignmentnum==0 and recordnum==0:
#         print(hsp.identities)

#print(len(readorthologs))



#for fasta in pangenome_sequences:
     #print((fasta.description.split(" "))[1])
     #if ((fasta.description.split(" "))[1]) in readorthologs:
    #     print((fasta.description.split(" "))[1])
         #print("So true!")

# OK so for some reason all these orthologs with
