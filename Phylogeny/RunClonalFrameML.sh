#!/bin/bash
# Use the concatenated core gene alignments output by CreateXMFA.sh 
# and ClonalFrameML to correct branch lengths for recombination 

source /home/acampbe/software/miniconda3/bin/activate TreeEnv
mkdir -p /home/acampbe/DFU/data/WGS_2020/Trees/clonalframeoutput
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/Trees/NoEpi100/RAxML_bestTree.RaxMLTreeNoRefs.newick /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/Trees/clonalframeoutput/clonalframe_tree.newick -xmfa_file true 
