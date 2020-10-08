
# Amy Campbell
# October 2020
# Identifying spatial organization of crtOPQMN genes in the DFU staph assemblies
# Do so via gff files

# All installed in xanthinEnv (Python3)
import pandas
import numpy
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('--Genomes_Exclude', default= "DORN1340,DORN1473,DORN672,DORN691,DORN701")
parser.add_argument('--geneNames', default = "group_2428,crtP,crtQ,crtM,crtN")
parser.add_argument('--gffDirPath', default="/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning/RoaryOutput/RoaryResults/gff_files")
parser.add_argument('--PresenceAbsencePath', default='/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning/RoaryOutput/RoaryResults/roaryoutput/gene_presence_absence.csv')
parser.add_argument('--OutputFilename', default='XanthinContigInfo')

args = parser.parse_args()

def IDstring(descript):
    #print(descript)
    splitlist = descript.split(';')
    return(str(splitlist[0]).replace("ID=", ""))

def CheckSameContig():

    return()

# List of assemblies to exclude (because we know them to have fragmented crtN genes)
genomes_exclude = (args.Genomes_Exclude).split(",")
# List of gff files not to look in
genomefiles_exclude = [g + ".gff" for g in genomes_exclude]

# In real life
gffdirpath = args.gffDirPath

#gff = pandas.read_csv(gffpath, sep='\t', header=0)
genenames = (args.geneNames).split(",")

pres_abs = pandas.read_csv(args.PresenceAbsencePath)
pres_abs = pres_abs.drop(genomes_exclude, axis=1)
pres_abs = pres_abs[pres_abs.Gene.isin(genenames)]

gffpathlist = os.listdir(gffdirpath)
gffpathlist = [g for g in gffpathlist if g not in genomefiles_exclude]
gffpathlist = [g for g in gffpathlist if g.endswith(".gff")]
genome_locDict = dict()

for gpath in gffpathlist:
    filepath= os.path.join(gffdirpath, gpath)
    DORN = gpath.replace(".gff","")

    ID_list = list(pres_abs[DORN])
    # Read in gff file (which is kind of a mess for any parser)
    gff = open(filepath, "r")

    # open a temporary file to write just the non-fasta business to
    tempfile = open("temporary.txt", "w")
    gffstring = gff.read()
    gffstring = gffstring.split("##FASTA")[0]
    tempfile.write(gffstring)
    tempfile.close()

    # Read it back in without the sequence stats or the fasta sequence (Via the temp file)
    pangff = pandas.read_csv("temporary.txt", sep='\t', header=None, comment="#", usecols=range(9))
    pangff.columns = ['contig', 'source','feature', 'start', 'end', 'score', 'strand', 'frame', 'description']

    # Format its ID column so we can find lines by ID
    pangff['ID'] = pangff.description.apply(lambda x: IDstring(x))

    generows = pangff[pangff.ID.isin(ID_list)]
    if generows['contig'].nunique() >> 1:
        genome_locDict[DORN] = ['False', generows['contig'].min(), generows['start'].min(), intgenerows['end'].max()]
    else:
        genome_locDict[DORN] = ['True', generows['contig'].min(), generows['start'].min(), generows['end'].max()]


FinalDF = pandas.DataFrame.from_dict(genome_locDict, orient='index')
FinalDF.columns = ['SameContig', 'Contig', 'StartPosition', 'EndPosition']
FinalDF.to_csv("XanthinContigInfo.csv", index=True)
