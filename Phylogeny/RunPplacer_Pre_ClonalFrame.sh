#!/bin/bash
# A posteriori placement of reference SA and S epidermidis onto the ML tree output by RaxML
# This tree is just for the purposes of presnations in march 21 (its preliminary, since it doesnt include clonalframeML step)

source /home/acampbe/software/miniconda3/bin/activate TreeEnv

#cp /home/acampbe/DFU/data/WGS_2020/Trees/NoRefs20/RAxML_info.RaxMLTreeNoRefs.newick /home/acampbe/DFU/data/WGS_2020/Trees/NoRefs20/RAxML_info.RaxMLTreeNoRefs.info

# Reference_Tree_File=/home/acampbe/DFU/data/WGS_2020/Trees/NoRefs20/RAxML_bestTree.RaxMLTreeNoRefs.newick

# remove pesky line from raxml 8.12 parameter output after renaming to RAxML_info.txt
# cp /home/acampbe/DFU/data/WGS_2020/Trees/NoRefs20/RAxML_info.RaxMLTreeNoRefs.newick /home/acampbe/DFU/data/WGS_2020/Trees/NoRefs20/RAxML_info.txt
# grep -v '^Partition: 0 with name:' RAxML_info.txt > RAxML_info.txt

TreeParams="/home/acampbe/DFU/data/WGS_2020/Trees/NoRefs20/RAxML_info.txt"

AllSeqs_Alignment="/home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/core_gene_alignment.fasta"

Reference_Tree_File=/home/acampbe/DFU/data/WGS_2020/Trees/NoRefs20/RAxML_bestTree.RaxMLTreeNoRefs.newick

# Run pplacer to make a mapping of these references a posteriori onto the RaxML generated tree
pplacer --mmap-file temp.txt -s $TreeParams -t $Reference_Tree_File $AllSeqs_Alignment

# Make a 'newick' file out of the mapping produced by pplacer
guppy tog -o core_gene_alignment.newick core_gene_alignment.jplace
mv core_gene_alignment.newick PPlacerTree207Isolates.newick
