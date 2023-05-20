#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforge/bin/activate ~/mambaforge/envs/pankmer

kmc_tools complex CC97_ingroups.txt
kmc_tools complex CC97_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC97_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC97_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC97_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC97_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC97_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC97_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC97_MarkerKmers.txt
