# Amy Campbell
# 12/2020

# Takes in the de novo assembled contigs
# Sticks all the circular fastas into one multifasta (where each sequence is delimited by DORNID_<contig #> s

import pandas
import os

inputfolder = "/home/acampbe/FinalContigs/"
outputfile = "/home/acampbe/DFU/data/WGS_2020/Plasmids/Plasmid_Multifasta.fasta"

filenames = os.listdir(inputfolder)
#os.mkdir(outputfolder)
output_table = open("/home/acampbe/DFU/data/WGS_2020/Plasmids/DORN_ShortRead_PlasmidSizes.txt", "w")

for filename in filenames:
    fasta = open(os.path.join(inputfolder, str(filename)), 'r')
    fastastring = fasta.read()
    stringlist = fastastring.split('>')
    plasmidlist = list(filter(lambda s: "circular=true" in s, stringlist))
    if len(plasmidlist) > 0:
        for p in plasmidlist:
            outputname = filename.replace("cleaned.fasta", "")
            length_p = ((p.split(" "))[1].split("length="))[1]

            output_table.write(str(filename.replace("_cleaned.fasta", str("_" + str(p_num))) + "," + length_p + "\n"))
            outputfile.write(str(">" + outputname + "_" + p + "\n"))
        outputfile.close()
