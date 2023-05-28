#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforgeOLD/bin/activate ~/mambaforgeOLD/envs/pankmer

kmc_tools complex CC72_ingroups.txt
kmc_tools complex CC72_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC72_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC72_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC72_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC72_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC72_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC72_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC72_MarkerKmers.txt
