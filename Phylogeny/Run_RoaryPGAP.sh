#!/bin/bash
# Updated Feb 2022

# bsub -e roarypgap.e -o roarypgap.o -n 8 -M 10240 -R "rusage [mem=10240] span[hosts=1]" sh Run_RoaryPGAP.sh 

source /home/acampbe/software/miniconda3/bin/activate prokenv
roary -e -z -p -8 -f /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022 /home/acampbe/DFU/data/WGS_2020/ReformattedPGAPNoSepi/*


