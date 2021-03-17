
from Bio import(SeqIO)
import pandas
import os

genepresence_table  = "/home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_bigMemory/gene_presence_absence.csv"
gene_alignments_folder = "/home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_bigMemory/pan_genome_sequences/"
outputfilepath = "/home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_bigMemory/core_genes.xmfa"
coregenelist = "/home/acampbe/DFUStrainsWGS/Phylogeny/core_gene_filelist.txt"


print(filelist)
genepresence_table = pandas.read_csv(genepresence_table)
genepresence_table = genepresence_table[genepresence_table.columns[14:len(genepresence_table.columns)]]
transposed_gene_presence = (genepresence_table.transpose())

output = open(outputfilepath, "w")
filelist = open(coregenelist, "r")
filelist = filelist.readlines()
filelist  = [fname.replace("\n", "") for fname in filelist]

i=0
for f in filelist:
    if not f.startswith('.'):
        multipath = os.path.join(gene_alignments_folder, f)
        genename = f.replace(".fa.aln", "")
        with open(multipath) as multifast:
            output.write("#" + str(genename) + "\n")
            for record in SeqIO.parse(multifast, "fasta"):
                row = transposed_gene_presence[transposed_gene_presence.apply(lambda r: r.str.contains(record.id, case=False).any(), axis=1)]
                genome = row.index.values[0]
                output.write(">" + str(genome) + "\n")
                output.write(str(record.seq) + "\n")
        output.write("=\n")
        print(genename)
        print(i)
        i=i+1

#genepresence_table[

#str.contains("APEBBLDI_01833").any()
#genome = [column for column in genepresence_table.columns if column.contains("APEBBLDI_01833").any()]

#print(genome)
