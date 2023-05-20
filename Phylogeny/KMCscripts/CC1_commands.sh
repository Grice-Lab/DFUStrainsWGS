#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforge/bin/activate ~/mambaforge/envs/pankmer

kmc_tools complex CC1_ingroups.txt
kmc_tools complex CC1_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC1_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC1_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC1_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC1_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC1_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC1_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC1_MarkerKmers.txt
