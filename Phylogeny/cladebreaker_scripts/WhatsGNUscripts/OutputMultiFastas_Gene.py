# Amy Campbell
# March 2023
# Preparing to run Mafft on the selected 'differentiating genes'
# Making a reference version of the gene (<whatever the first DORN it encounters is>_genename.fasta)
# Then making a multifasta of everyone else's version of the gene in the format, for example of Patient176:
# Patient176_metI.fasta:
# > DORN1646_metI
# ACTGAGCT....
# > GCA_000555645.1_metI
# ACTGAGCT....

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas
import os
import sys

# Takes in:
############

# - Path to the FFNs for input (arg 2, should be folder containing only the .ffns you wanna use of the format <genomename>.ffn)
    # /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/FFNs
# - Path to the roary CSV path (arg 3)
    # /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/roary_alignment/results/gene_presence_absence.csv
# - Path to unique genes list output by UniqueWithinGenome.py (arg 4)
    # /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/Patient176_unique.txt
# - Folder to put the mafft calls script
    # /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/



# puts out "GeneFastas" folder containing fastas with each gene's sequence in each genome
# as well as one randomly chosen DORN reference genome's
# as well as a script called MafftAlignGenes.sh
#########################################################################################

argsinput = sys.argv
if len(argsinput) <5:
    print("Not enough input files actually")
    exit()
else:
    FFNpath = argsinput[1]
    roaryCSVpath = argsinput[2]
    UniqueGenesList = argsinput[3]
    OutputScriptDir = argsinput[4]

# Read in the reference genome's ffn (nt sequences for each )

OutputScript=os.path.join(OutputScriptDir, "MafftAlignGenes.sh")

IncludedGenomes = (os.listdir(FFNpath))

IncludedDORNs = list(filter(lambda i: 'DORN' in i, IncludedGenomes))

IncludedGenomeNames = list(map(lambda x: x.replace(".ffn", ""),IncludedGenomes))

ReferenceFFN = IncludedDORNs[0]
ReferenceName=ReferenceFFN.replace(".ffn", "")
ReferencePath=os.path.join(FFNpath, ReferenceFFN)

# Make a big ugly dictionary to pull from
AllDicts = dict()
for G in IncludedGenomes:
    path_fasta=os.path.join(FFNpath, G)
    GenomeString=G.replace(".ffn", "")
    AllDicts[GenomeString] = SeqIO.to_dict(SeqIO.parse(path_fasta, "fasta"))

# Read in the Unique Genes' names
ListOrthologs = open(UniqueGenesList, "r").read().split('\n')
ListOrthologs.remove('')

# Read in the gene presence absence csv from roary as a pandas df
OrthDF = pandas.read_csv(roaryCSVpath)
OrthDF = OrthDF[OrthDF['Gene'].isin(ListOrthologs)]
#SeqIO.write(sequences,handle,"fasta")

# make an output directory
OutputKey = os.path.dirname(UniqueGenesList)
OutputDir = os.path.join(OutputKey, "GeneFastas")
os.mkdir(OutputDir)

ReferenceOutputPathPrefix = os.path.join(OutputDir, str(os.path.basename(OutputKey)) + ReferenceName )
MultiOutputPathPrefix = os.path.join(OutputDir, str(os.path.basename(OutputKey))  )

# example of making output file
#print(str(ReferenceOutputPathPrefix)+ "_something.fasta")
#with open("example.fasta", "w") as output_handle:
#    SeqIO.write(sequences, output_handle, "fasta")


FinalListOrths = []
for o in ListOrthologs:
    rowItem = (OrthDF[OrthDF.Gene==str(o)])

    ReferenceProtein = str(rowItem[ReferenceName].item())
    if ReferenceProtein=='nan':
        print(str(o) + " is missing from reference genome " + str(ReferenceName) + "! Removing from Orthologs List.")
    else:
        dictionary_item = (AllDicts[ReferenceName][ReferenceProtein]) # [str(ReferenceProtein)].seq )
        with open(str(ReferenceOutputPathPrefix) + "_" +  str(o) + "_" + ReferenceName + ".fasta", "w") as refotpt:
            dictionary_item.id = ReferenceName
            SeqIO.write(dictionary_item,refotpt, "fasta")

        with open(str(MultiOutputPathPrefix) + "_" + str(o) + ".fasta", "w") as multioutput:
            for g in IncludedGenomeNames:

                ProteinID = str(rowItem[g].item())
                # If any one genome is *missing* this gene (weird but possible I guess), remove it from teh list and say so.
                if ProteinID == 'nan':
                    print(str(o) + " is missing from " + str(g) + "! Removing from Orthologs List.")
                else:
                    FinalListOrths.append(o)
                    dictionary_item = (AllDicts[g][ProteinID])
                    dictionary_item.id = g
                    SeqIO.write(dictionary_item, multioutput, "fasta")

ListOrths = open(str(MultiOutputPathPrefix + "_FinalOrthologs.txt"), "w")
for F in FinalListOrths:
    ListOrths.write(str(F) + "\n")

ListOrths.close()

with open(OutputScript, "w") as scriptoutput:
    scriptoutput.write("#!bin/bash\n")
    scriptoutput.write("# Mafft Alignment on Genes \n")
    scriptoutput.write("################\n")
    for f in FinalListOrths:
        outputstr = "mafft --localpair --addfragments" +  str(MultiOutputPathPrefix) + "_" + str(f) + ".fasta " + str(ReferenceOutputPathPrefix) + "_" +  str(f) + "_" + ReferenceName + ".fasta  > MafftOutput_" + str(f) + "\n"
        scriptoutput.write(outputstr)
