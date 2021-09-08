#!/bin/bash
# Align and compare whole genomes from patient 141

source /home/acampbe/software/miniconda3/bin/activate Patient141


# nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff DORN925_DORN1088Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
	/home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
	/home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_925_1000 DORN925_DORN1000Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_925_1037 DORN925_DORN1037Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_925_1194 DORN925_DORN1194Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_925_1061 DORN1925_DORN1061Comparison




nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1037_1000 DORN1037_DORN1000Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1061_1000 DORN1061_DORN1000Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1000_1088 DORN1000_DORN1088Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1194_1000 DORN1194_DORN1000Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1037_1061 DORN1037_DORN1061Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1037_1088 DORN1037_DORN1088Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1097_1194 DORN1037_DORN1194Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1061_1088 DORN1061_DORN1088Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1061_1194 DORN1061_DORN1194

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta \
	/home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_1194_1088 DORN1194_DORN1088
