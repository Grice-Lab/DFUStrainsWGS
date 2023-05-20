#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforge/bin/activate ~/mambaforge/envs/pankmer

kmc_tools complex CC12_ingroups.txt
kmc_tools complex CC12_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC12_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC12_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC12_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC12_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC12_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC12_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC12_MarkerKmers.txt
