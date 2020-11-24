
# Amy Campbell
# October 2020
# Use the loci identified by ID_Sxanthin_Loci.R to extract the crtOPQMN operon

import pandas
import os
inputfolder = "/project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Cleaned_Assemblies/FinalContigs/"

inputfolder="/home/acampbe/FinalContigs/"

filenames = ["DORN315_cleaned.fasta","DORN2144_cleaned.fasta","DORN2187_cleaned.fasta",
"DORN2083_cleaned.fasta","DORN2089_cleaned.fasta","DORN2123_cleaned.fasta",
"DORN1702_cleaned.fasta","DORN1779_cleaned.fasta","DORN1761_cleaned.fasta",
"DORN1740_cleaned.fasta","DORN1289_cleaned.fasta","DORN1285_cleaned.fasta",
"DORN1881_cleaned.fasta","DORN339_cleaned.fasta","DORN333_cleaned.fasta",
"DORN317_cleaned.fasta","DORN280_cleaned.fasta","DORN300_cleaned.fasta",
"DORN283_cleaned.fasta","DORN781_cleaned.fasta","DORN767_cleaned.fasta",
"DORN869_cleaned.fasta","DORN825_cleaned.fasta","DORN2148_cleaned.fasta",
"DORN2213_cleaned.fasta","DORN2179_cleaned.fasta","DORN2136_cleaned.fasta",
"DORN2105_cleaned.fasta","DORN2149_cleaned.fasta","DORN2150_cleaned.fasta",
"DORN2137_cleaned.fasta","DORN498_cleaned.fasta","DORN460_cleaned.fasta",
"DORN489_cleaned.fasta","DORN467_cleaned.fasta","DORN564_cleaned.fasta",
"DORN459_cleaned.fasta","DORN429_cleaned.fasta","DORN468_cleaned.fasta",
"DORN1465_cleaned.fasta","DORN1447_cleaned.fasta","DORN1455_cleaned.fasta",
"DORN2064_cleaned.fasta","DORN2059_cleaned.fasta","DORN1373_cleaned.fasta",
"DORN1368_cleaned.fasta","DORN1358_cleaned.fasta","DORN359_cleaned.fasta",
"DORN1943_cleaned.fasta","DORN377_cleaned.fasta","DORN383_cleaned.fasta",
"DORN1858_cleaned.fasta","DORN1885_cleaned.fasta","DORN1863_cleaned.fasta",
"DORN1834_cleaned.fasta","DORN1808_cleaned.fasta","DORN1923_cleaned.fasta",
"DORN1644_cleaned.fasta","DORN1695_cleaned.fasta","DORN1585_cleaned.fasta",
"DORN1663_cleaned.fasta","DORN1731_cleaned.fasta","DORN1643_cleaned.fasta",
"DORN1818_cleaned.fasta","DORN1819_cleaned.fasta","DORN1776_cleaned.fasta",
"DORN807_cleaned.fasta","DORN1844_cleaned.fasta","DORN2127_cleaned.fasta",
"DORN1869_cleaned.fasta","DORN2034_cleaned.fasta","DORN2139_cleaned.fasta",
"DORN2205_cleaned.fasta","DORN2166_cleaned.fasta","DORN1849_cleaned.fasta",
"DORN1902_cleaned.fasta","DORN1952_cleaned.fasta","DORN2004_cleaned.fasta",
"DORN1968_cleaned.fasta","DORN2075_cleaned.fasta","DORN1176_cleaned.fasta",
"DORN1829_cleaned.fasta"]


output_table = open("PlasmidSizes.txt", "w")
for filename in filenames:
    fasta = open(os.path.join(inputfolder, str(filename)), 'r')
    fastastring = fasta.read()
    stringlist = fastastring.split('>')
    plasmidlist = list(filter(lambda s: "circular=true" in s, stringlist))
    if len(plasmidlist) > 0:
        outputname = filename.replace("cleaned.fasta", "circular.fasta")
        outputfile = open(str("/home/acampbe/CC5Plasmids/" + outputname), "w")


        for p in plasmidlist:
            length_p = ((p.split(" "))[1].split("length="))[1]
            output_table.write(str(filename.replace("_cleaned.fasta", "")+ "," + length_p + "\n"))
            outputfile.write(str(">" + p + "\n"))

        outputfile.close()
