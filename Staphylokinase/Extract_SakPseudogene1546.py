# Amy Campbell
# Extract what I think is a pseudogene of sak from DORN1546 

from Bio import SeqIO


# e.g. /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta
DORN1546Path = "/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1546_Final.fasta"

# e.g. /home/acampbe/DFU/data/WGS_2020/ONP/Patient141/DORN925Phage.fasta
outputPath = "sakPseudogeneDORN1546.fasta"

contigs = list(SeqIO.parse(DORN1546Path, "fasta"))
contig1 = str(contigs[0].seq)

geneseq = contig1[1712989:1714083]


otpt = open(outputPath, "w")
otpt.write(">sakPseudogene")
otpt.write(geneseq)
otpt.close()
