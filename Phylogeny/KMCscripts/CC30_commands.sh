#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforgeOLD/bin/activate ~/mambaforgeOLD/envs/pankmer

kmc_tools complex CC30_ingroups.txt
kmc_tools complex CC30_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC30_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC30_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC30_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC30_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC30_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC30_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC30_MarkerKmers.txt
