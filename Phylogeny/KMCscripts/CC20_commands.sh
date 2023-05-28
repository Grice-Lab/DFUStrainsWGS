#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforgeOLD/bin/activate ~/mambaforgeOLD/envs/pankmer

kmc_tools complex CC20_ingroups.txt
kmc_tools complex CC20_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC20_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC20_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC20_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC20_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC20_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC20_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC20_MarkerKmers.txt
