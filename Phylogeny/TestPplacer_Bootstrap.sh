#!/bin/bash
# A posteriori placement of reference SA and S epidermidis onto the ML tree output by RaxML

source /home/acampbe/software/miniconda3/bin/activate TreeEnv

# RaxML needs the .fasta extension for the alignment file 
# cp /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/core_gene_alignment.aln /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/core_gene_alignment.fasta
mv /home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_info.RaxMLTreeNoRefs.newick /home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_info.NoTreeRefs.info

TreeParams=/home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_info.NoTreeRefs.info
#RAxML_info.RaxMLTreeNoRefs.newick

#RAxML_info.NoTreeRefs.info
Reference_Tree_File=/home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_bestTree.RaxMLTreeNoRefs.newick

#TreeParams=/home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_info.RaxMLTreeNoRefs.info

AllSeqs_Alignment=/home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/core_gene_alignment.fasta

# Run pplacer to make a mapping of these references a posteriori onto the RaxML generated tree
pplacer --mmap-file temp.txt -s $TreeParams -t $Reference_Tree_File $AllSeqs_Alignment

# Make a 'newick' file out of the mapping produced by pplacer
guppy tog -o core_gene_alignment.newick core_gene_alignment.jplace

bootstrapLabeled="/home/acampbe/DFU/data/WGS_2020/Trees/PartitionedNoRefs/RAxML_bipartitionsBranchLabels.RaxMLTreeNoRefsParitioned"

