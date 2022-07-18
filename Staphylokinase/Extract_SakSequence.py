# Amy Campbell
# Extract sak gene sequence for Amelia to check her PCR primers

from Bio import SeqIO


# e.g. /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta
DORN925Path = "/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta"

# e.g. /home/acampbe/DFU/data/WGS_2020/ONP/Patient141/DORN925sak.fasta
outputPath = "DORN925sak.fasta"

contigs = list(SeqIO.parse(DORN925Path, "fasta"))
contig1 = str(contigs[0].seq)

geneseq = contig1[22091:22579]


otpt = open(outputPath, "w")
otpt.write(">sak925")
otpt.write(geneseq)
otpt.close()
