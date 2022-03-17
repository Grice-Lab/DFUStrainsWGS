#!/bin/bash
# Align and compare whole genomes from patient 141
# specifically in this case compare each genome to 925

source /home/acampbe/software/miniconda3/bin/activate Patient141


nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN929_Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_929 DORN925_DORN929Comparison

