
from Bio import SeqIO
import sys


# Marker 1
#  1554050:1554550 (SNP 4)

# Marker 2
#  2352094:2352594 (SNP 16)

# Marker 3
#  2729951:2730451 (SNP 13)

# Marker 4
#  2752692:2753200 (SNP 10)

# e.g. /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta
DORN925Path = sys.argv[1]

# e.g. /home/acampbe/DFU/data/WGS_2020/ONP/Patient141/Markers/XanthinMarkers.fasta
outputPath = sys.argv[2]

contigs = list(SeqIO.parse(DORN925Path, "fasta"))
contig1 = str(contigs[0].seq)

marker1 = contig1[1554050:1554550]

marker2 = contig1[2352094:2352594]

marker3 = contig1[2729951:2730451]

marker4 = contig1[2752692:2753200]

otpt = open(outputPath, "w")
otpt.write(">Marker1\n")
otpt.write(marker1)
otpt.write("\n>Marker2\n")
otpt.write(marker2)
otpt.write("\n>Marker3\n")
otpt.write(marker3)
otpt.write("\n>Marker4\n")
otpt.write(marker4)

otpt.close()
