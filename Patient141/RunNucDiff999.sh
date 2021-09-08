#!/bin/bash
# Align and compare whole genomes from patient 141

source /home/acampbe/software/miniconda3/bin/activate Patient141


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_925_999 DORN925_DORN999Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_925_976 DORN925_DORN976Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_976_999 DORN976_DORN999Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_976_1000 DORN976_DORN1000Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_976_1037 DORN976_DORN1037Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_976_1061 DORN976_DORN1061Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_976_1088 DORN976_DORN1088Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_976_1194 DORN976_DORN1194Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_999_1000 DORN999_DORN1000Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_999_1037 DORN999_DORN1037Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_999_1061 DORN999_DORN1061Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_999_1194 DORN999_DORN1194Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/NucDiff_999_1088 DORN999_DORN1088Comparison
