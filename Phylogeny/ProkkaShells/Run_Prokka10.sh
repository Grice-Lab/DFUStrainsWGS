#!/bin/bash
mkdir -p /home/acampbe/DFU/data/WGS_2020/RoaryResults
source /home/acampbe/software/miniconda3/bin/activate prokenv
prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/DORN1285 --force --prefix DORN1285 --genus Staphylococcus /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Cleaned_Assemblies/FinalContigs/DORN1285_cleaned.fasta
prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/DORN1289 --force --prefix DORN1289 --genus Staphylococcus /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Cleaned_Assemblies/FinalContigs/DORN1289_cleaned.fasta
prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/DORN1844 --force --prefix DORN1844 --genus Staphylococcus /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Cleaned_Assemblies/FinalContigs/DORN1844_cleaned.fasta
