
from Bio import SeqIO
import sys


# Marker 1
#  1554150:1554450 (SNP 4)

# Marker 2
#  2352194:2352494 (SNP 16)

# Marker 3
#  2730051:2730351 (SNP 13)

# Marker 4
#  2752792:2753100 (SNP 10)

# e.g. /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta
DORN925Path = sys.argv[1]

# e.g. /home/acampbe/DFU/data/WGS_2020/ONP/Patient141/Markers/XanthinMarkers.fasta
outputPath = sys.argv[2]

contigs = list(SeqIO.parse(DORN925Path, "fasta"))
contig1 = str(contigs[0].seq)

marker1 = contig1[1554150:1554450]

marker2 = contig1[2352194:2352494]

marker3 = contig1[2730051:2730351]

marker4 = contig1[2752792:2753100]

otpt = open(outputPath, "w")
otpt.write("> Marker1\n")
otpt.write(marker1)
otpt.write("\n> Marker2\n")
otpt.write(marker2)
otpt.write("\n> Marker3\n")
otpt.write(marker3)
otpt.write("\n> Marker4\n")
otpt.write(marker4)

otpt.close()
