# Amy Campbell Feb 2022 
#
# Modified from MakeRoaryShells.py from Feb 2020 
# Makes shells to run prokka on chunks of the DORN fastas at once
# takes in the Genus name (default S. aureus since I'm using this for the DORNs)
# as well as info about the input directory, output directory for the shellscripts
# to put vs. the data to put out (like the .gff files produced by prokka)

# example call


# inputdir="/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/"
# outdirshells="/home/acampbe/DFUStrainsWGS/Phylogeny/ProkkaShells/"
# python3 MakeProkkaShells.py --inputDirPath $inputdir --outputdirshells $outdirshells --outputdir_prokka /home/acampbe/DFU/data/WGS_2020/completed_prokka --conda_activatepath /home/acampbe/software/miniconda3/bin/activate --InputExtension "_Final.fasta" --nshells 10

import argparse
import os

# Process user input
parser = argparse.ArgumentParser(description = "Process user input")

parser.add_argument('--inputDirPath', type=str)

parser.add_argument('--outputdirshells', type=str)

parser.add_argument('--outputdir_prokka', type=str)

parser.add_argument('--conda_activatepath',
                    type=str, help="The path at which your /bin/conda/activate is located")

parser.add_argument('--nshells', type=int, help="The number of shell scripts you'd like to divide your calls to Prokka into." +
                                                 " If  you only have 10 scripts, say, you could enter 1 for this. But if you have 300...")
parser.add_argument('--InputExtension',default="_contigs.fasta", type=str, help="This script assumes the fastq filenames are in the format of <sampleID>_contigs.fasta")

parser.add_argument('--shellscriptname', default="Run_Prokka", type=str,
help="By default, shellscript filenames are in the format of Run_Prokka_<#>.sh." + "Use this flag to change Run_Prokka to something else")

parser.add_argument('--genusname', default="Staphylococcus", type=str)

parser.add_argument('--commonString', default="DORN", type=str, help="Shared string between all filenames.")

args = parser.parse_args()

# Save args as variables
#########################
contigspath = args.inputDirPath
outputdirprokka = str(args.outputdir_prokka)
outputdirshells= str(args.outputdirshells)
numshells = int(args.nshells)
common = str(args.commonString)
inputextension = str(args.InputExtension)
condapath = str(args.conda_activatepath)
shellname = str(args.shellscriptname)
genus = str(args.genusname)

# Get list of fastas
fnamelist = os.listdir(contigspath)
fnamelist = list(filter(lambda k: inputextension in k, fnamelist))
fnamelist = list(filter(lambda k: common in k, fnamelist))


# Lazily distribute the samples across shell scripts

subsize = len(fnamelist) // numshells

if numshells==1:
    chunk_array = [fnamelist]
else:
    chunk_array=[]
    starting = 0
    for i in range(numshells):
        chunk_array.append(fnamelist[starting:(starting+subsize)])
        starting=starting+subsize
    chunk_array.append(fnamelist[starting:len(fnamelist)])

for c in range(len(chunk_array)):
    outputfile = os.path.join(outputdirshells, str(shellname + str(c) + ".sh"))
    otpt = open(outputfile, "w")
    otpt.write("#!/bin/bash\n")
    otpt.write(str("mkdir -p "+ str(outputdirprokka)+ "\n"))
    otpt.write(str("source " + condapath + " prokenv\n"))

    # For each assembly in this 'chunk', add to the shellscript it's in,
    # Make a new file for it.
    output_gff_folder = os.path.join(str(outputdirprokka), "gff_files")
    for j in chunk_array[c]:
        genomestring = j.replace(inputextension, "")
        inputstring = os.path.join(contigspath, j)
        otpt.write(str("prokka --outdir " + str(output_gff_folder)+"/"+ genomestring+ " --force "+ "--prefix "+ str(genomestring) + " --genus " + genus + " " + inputstring + "\n"))
    otpt.close()














#
