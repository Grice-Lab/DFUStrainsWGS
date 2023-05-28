#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforgeOLD/bin/activate ~/mambaforgeOLD/envs/pankmer

kmc_tools complex CC59_ingroups.txt
kmc_tools complex CC59_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC59_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC59_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC59_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC59_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC59_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC59_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC59_MarkerKmers.txt
