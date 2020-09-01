#!/bin/bash
# Runs RaxML max likelihood tree construction
#  on the core genome alignment output by Roary (this is really a rough estimate of this alignment as compared to Mugsy output)
source /home/acampbe/software/miniconda3/bin/activate TreeEnv

# Run July 26 -- RaxML on core gene alignment by Roary with PRANK (no reference genomes; we'll add those a posteirori)
# Estimate tree
raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln -n RaxMLTreeNoRefs.newick 


# THIS NEXT PART IS ACTUALLY UNNECESSARY 
# Outuput params of that tree
# raxmlHPC -f e -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln -t RAxML_bestTree.RaxMLTreeNoRefs.newick -n PARAMS



