
# Amy Campbell
# Nov 2020
# extract and save sizes of plasmids circularized by plasmidspades for CC8 complex

import pandas
import os
inputfolder = "/project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Cleaned_Assemblies/FinalContigs/"

inputfolder="/home/acampbe/FinalContigs/"

filenames = ["DORN1755_cleaned.fasta",
"DORN1691_cleaned.fasta",
"DORN1657_cleaned.fasta",
"DORN1661_cleaned.fasta",
"DORN1686_cleaned.fasta",
"DORN434_cleaned.fasta",
"DORN1207_cleaned.fasta",
"DORN1253_cleaned.fasta",
"DORN1219_cleaned.fasta",
"DORN657_cleaned.fasta",
"DORN717_cleaned.fasta",
"DORN76_cleaned.fasta",
"DORN62_cleaned.fasta",
"DORN105_cleaned.fasta",
"DORN47_cleaned.fasta",
"DORN56_cleaned.fasta",
"DORN1502_cleaned.fasta",
"DORN1499_cleaned.fasta",
"DORN1743_cleaned.fasta",
"DORN1782_cleaned.fasta",
"DORN1747_cleaned.fasta",
"DORN1729_cleaned.fasta",
"DORN1765_cleaned.fasta",
"DORN1214_cleaned.fasta",
"DORN2221_cleaned.fasta",
"DORN2058_cleaned.fasta",
"DORN2093_cleaned.fasta",
"DORN2130_cleaned.fasta",
"DORN1946_cleaned.fasta",
"DORN1961_cleaned.fasta",
"DORN848_cleaned.fasta",
"DORN892_cleaned.fasta",
"DORN788_cleaned.fasta",
"DORN731_cleaned.fasta",
"DORN936_cleaned.fasta",
"DORN1006_cleaned.fasta",
"DORN973_cleaned.fasta"]

output_table = open("CC8_PlasmidSizes.txt", "w")
for filename in filenames:
    fasta = open(os.path.join(inputfolder, str(filename)), 'r')
    fastastring = fasta.read()
    stringlist = fastastring.split('>')
    plasmidlist = list(filter(lambda s: "circular=true" in s, stringlist))
    if len(plasmidlist) > 0:
        outputname = filename.replace("cleaned.fasta", "circular.fasta")
        outputfile = open(str("/home/acampbe/CC8Plasmids/" + outputname), "w")

        for p in plasmidlist:
            length_p = ((p.split(" "))[1].split("length="))[1]
            output_table.write(str(filename.replace("_cleaned.fasta", "")+ "," + length_p + "\n"))
            outputfile.write(str(">" + p + "\n"))

        outputfile.close()
