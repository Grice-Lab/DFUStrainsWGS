#!/bin/bash
# Align and compare whole genomes of DORN925 vs DORN1088

source /home/acampbe/software/miniconda3/bin/activate Patient141

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff DORN925_DORN1088Comparison

