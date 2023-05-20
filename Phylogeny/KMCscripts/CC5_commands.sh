#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforge/bin/activate ~/mambaforge/envs/pankmer

kmc_tools complex CC5_ingroups.txt
kmc_tools complex CC5_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC5_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC5_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC5_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC5_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC5_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC5_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC5_MarkerKmers.txt
