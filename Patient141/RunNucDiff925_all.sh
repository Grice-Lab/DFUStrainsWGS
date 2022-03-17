#!/bin/bash
# Align and compare whole genomes from patient 141
# specifically in this case compare each genome to 925

source /home/acampbe/software/miniconda3/bin/activate Patient141

mkdir -p /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1000.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1000 DORN925_DORN1000Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN976.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_976 DORN925_DORN976Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN999.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_999 DORN925_DORN999Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1194.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1194 DORN925_DORN1194Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1061.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1061 DORN925_DORN1061Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1088.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1088 DORN925_DORN1088Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1038.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1038 DORN925_DORN1038Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN881.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_881 DORN925_DORN881Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN880.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_880 DORN925_DORN880Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN1037.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1037 DORN925_DORN1037Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN933.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_933 DORN925_DORN933Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN925.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/HybridAssemblies/DORN882.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_882 DORN925_DORN882Comparison


