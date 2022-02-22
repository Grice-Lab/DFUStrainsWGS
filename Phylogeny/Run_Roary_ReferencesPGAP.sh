#!/bin/bash
# bsub -e roarySepiPGAP.e -o roarySepiPGAP.o  -n 4 -M 10240 -R "rusage [mem=10240] span[hosts=1]" sh Run_Roary_ReferencesPGAP.sh

source /home/acampbe/software/miniconda3/bin/activate prokenv
roary -e -p -4 -f /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022_Epi /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/*

