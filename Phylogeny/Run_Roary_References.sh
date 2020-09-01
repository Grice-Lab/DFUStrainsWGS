#!/bin/bash
source /home/acampbe/software/miniconda3/bin/activate prokenv
roary -e -p -8 -f /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs /home/acampbe/DFU/data/WGS_2020/RoaryResults/gffs_wRefs/*
cp /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/core_gene_alignment.aln /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/core_gene_alignment.fasta
