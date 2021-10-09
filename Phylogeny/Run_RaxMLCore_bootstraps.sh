#!/bin/bash
# Runs RaxML max likelihood tree construction with bootstrapping

source /home/acampbe/software/miniconda3/bin/activate TreeEnv

mkdir -p /home/acampbe/DFU/data/WGS_2020/Trees/NoRefsBootstrapped 


# Run Feb 2021 -- RaxML on core gene alignment by Roary with PRANK (no reference genomes; we'll add those a posteirori) WITH bootstrapping
# Estimate tree
# raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln -n RaxMLTreeNoRefs.newick 

#raxmlHPC -m GTRGAMMA -p 19104 -x 19104 -# 100 -s /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln -n RaxMLTreeNoRefsBS.newick
#mv *.newick /home/acampbe/DFU/data/WGS_2020/Trees/NoRefsBootstrapped/

raxmlHPC -m GTRGAMMA -p 19104 -x 19104 -# 10 -s /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln -n RaxMLTreeNoRefsTEST.newick
