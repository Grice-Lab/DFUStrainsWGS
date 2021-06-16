#!/bin/bash
# Amy Campbell
# 03/2021
# take gene-by-gene alignment outputs by prank with roary call including -z (keep intermediate files)
# concatenate them for input into clonalframeML

#source /home/acampbe/software/miniconda3/bin/activate python3env

source /home/acampbe/software/miniconda3/bin/activate R_envir

# Make core gene list
Rscript List_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/ core_gene_filelist.txt 231

source /home/acampbe/software/miniconda3/bin/activate python3env

# Concatenate them
python3 Create_XMFA_File.py
