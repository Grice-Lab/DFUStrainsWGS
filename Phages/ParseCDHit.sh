#!/bin/bash
# Amy Campbell
# Parse the output of CD-hit-est
# to make a new multifasta with nicer names (Phage1...Phage68)
# and also keep track of which genomes had or didn't have that phage 
# cluster (complete according to CheckV), prior to making a mapping-based
# estimate

source ~/mambaforge/bin/activate DBSCANSWA

clstr_file=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/ClusteredPhages90.clstr
fastainput=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/ClusteredPhages90
GenomesFolder=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates
CSVoutput=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/PhagePresenceAbsence.csv
fastaoutput=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/ClusteredPhagesRenamed.fasta

python3 ParseCDhit.py $clstr_file $fastainput $GenomesFolder $CSVoutput $fastaoutput
