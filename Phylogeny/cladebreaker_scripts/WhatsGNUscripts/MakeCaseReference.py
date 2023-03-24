# Amy Campbell
# March 2023
# Makes a reference fasta of marker genes to align metagenomes against
# Takes in:
# 1. Path to 'GeneFastas' which includes the reference genome's version of each gene
#   (e.g. /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/GeneFastas)
#   which contains Patient176DORN1863_hutG_DORN1863.fasta and others
# 2. Path to the list of genes to actually use (Patient176_Variants_FinalGeneList.txt for example)
# 3. Patient ID used in fasta names (e.g. "Patient176")
# 4. Output file path

# reads in the list of genes to use

# For each gene,
# Reads in the FASTA <PatientName(e.g. Patient176)><referenceName>_<genename>_<referenceName>.fasta
# as a SeqIO record
# Updates the id to be the gene name
# writes it out to the output file
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

# Example call:
# GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/GeneFastas"
# FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/Patient176_Variants_FinalGeneList.txt"
# PatientID="Patient176"
# OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/"
# ReferenceGenome="DORN1863"
# python3 MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

argsinput = sys.argv
if len(argsinput) <6:
    print("Not enough input files actually")
    exit()
else:
    GeneFastasDir = argsinput[1]
    FinalGeneList = argsinput[2]
    PatientID = argsinput[3]
    OutputfilePath = argsinput[4]
    RefGenomeName=argsinput[5]



# Open list of genes
IncludedGenes = open(FinalGeneList, "r").read().split('\n')
IncludedGenes.remove('')

OtptFasta = open(OutputfilePath, "w")
for g in IncludedGenes:
    pathToOpen = os.path.join(GeneFastasDir, str(PatientID) +RefGenomeName + "_"+g+ "_"+RefGenomeName+".fasta" )
    SeqObj = next(SeqIO.parse(pathToOpen, "fasta"))

    OutputObj=SeqIO.SeqRecord(Seq(SeqObj.seq), id=g, description=RefGenomeName)
    SeqIO.write(OutputObj, OtptFasta, "fasta")

OtptFasta.close()
