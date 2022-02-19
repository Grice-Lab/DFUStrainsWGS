#!/bin/bash
mkdir -p /home/acampbe/DFU/data/WGS_2020/completed_prokka
source /home/acampbe/software/miniconda3/bin/activate prokenv
prokka --outdir /home/acampbe/DFU/data/WGS_2020/completed_prokka/gff_files/DORN383 --force --prefix DORN383 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/PGAPformattedIsolates/DFUIsolates/DORN383_pgap.fasta
