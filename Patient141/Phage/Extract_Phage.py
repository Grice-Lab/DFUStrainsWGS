# Amy Campbell
# Extracting phage sequences from DORN925 for Patric annotation

from Bio import SeqIO
import sys

# python3 /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.f$

# e.g. /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta
DORN925Path = sys.argv[1]

# e.g. /home/acampbe/DFU/data/WGS_2020/ONP/Patient141/DORN925Phage.fasta
outputPath = sys.argv[2]

contigs = list(SeqIO.parse(DORN925Path, "fasta"))
contig1 = str(contigs[0].seq)

phageend = contig1[2789823:2806769]
phagebeginning=contig1[1:25619]

otpt = open(outputPath, "w")
otpt.write(">DORN925phage1\n")
otpt.write(phageend)
otpt.write("\n>DORN925phage2\n")
otpt.write(phagebeginning)
otpt.close()
