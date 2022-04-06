# Extract +/- 50 bp from either end of each SNP and variant 

source /home/acampbe/software/miniconda3/bin/activate prokenv

python3 ExtractMetaMarkers141_Extended.py /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta /home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/XanthinMarkersExtended.fasta

