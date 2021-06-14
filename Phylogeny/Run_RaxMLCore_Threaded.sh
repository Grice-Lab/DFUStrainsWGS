#!/bin/bash
# Runs RaxML max likelihood tree construction
# on the core genome alignment output by Roary via PRANK
# This script in particular uses 8 threads with raxmlHPC-PTHREADS
# so has to be run with 
# bsub -e runraxmlthreaded.e -o runraxmlthreaded.o -n 4 -M 10240 -R "rusage [mem=10240] span[hosts=1]" sh Run_RaxMLCore_Threaded.sh 

source /home/acampbe/software/miniconda3/bin/activate TreeEnv

mkdir -p /home/acampbe/DFU/data/WGS_2020/Trees/NoEpi100
mkdir -p /home/acampbe/DFU/data/WGS_2020/Trees/NoEpiBootstrapped1000/

# Estimate best scoring ML tree from 20 random starts 
raxmlHPC-PTHREADS -T 8 -m GTRGAMMA -p 19104 -# 100 -s /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln -n RaxMLTreeNoRefs.newick 

mv *newick* /home/acampbe/DFU/data/WGS_2020/Trees/NoEpi100/

# Perform 1000 bootstraps
raxmlHPC-PTHREADS -T 8 -m GTRGAMMA -p 19104 -x 19104 -# 1000 -s /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln -n RaxMLTreeNoRefsBS.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/Trees/NoEpiBootstrapped1000/

# Use the bootstrapped topologies to draw partitions on the best ML tree 
raxmlHPC -m GTRGAMMA -p 19104 -f b -t /home/acampbe/DFU/data/WGS_2020/Trees/NoEpi100/RAxML_bestTree.RaxMLTreeNoRefs.newick -z /home/acampbe/DFU/data/WGS_2020/Trees/NoEpiBootstrapped1000/RAxML_bootstrap.RaxMLTreeNoRefsBS.newick -n RaxMLTreeParitioned

