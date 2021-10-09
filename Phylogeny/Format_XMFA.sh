#!/bin/bash
# Amy Campbell
# 03/2021
# Format the gene-by-gene alignments output by the roary -z option into the 
# pan_genome_sequences/ folder into an .xmfa file in which multifasta blocks
# (entire file contents) are separated by a line with '=' (apparently the file
# also has to END in a '=' line; not sure why but it's some parsing thing)
# Note: while pan_genome_sequences/ has all the genes in the pan genome aligned, 
# I need to just grab the core genes that were used to build the core gene alignment (and therefore used to build the tree)

filelist="core_gene_filelist.txt"

fileLocation="/home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_bigMemory/pan_genome_sequences/"

outputfile="/home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_genes.xmfa"

# Make output file
cat > $outputfile 

# read filelist, adding the contents of each to the output file separated by "=" 

while IFS= read -r line;
do
	cat $fileLocation$line >> $outputfile
	echo "=" >> $outputfile

done < "$filelist"



