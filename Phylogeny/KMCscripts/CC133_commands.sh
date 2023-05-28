#!bin/bash
# Performing KMC operations for clade-specific markers

source ~/mambaforgeOLD/bin/activate ~/mambaforgeOLD/envs/pankmer

kmc_tools complex CC133_ingroups.txt
kmc_tools complex CC133_outgroups.txt
kmc_tools simple /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC133_intersect /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC133_Union kmers_subtract /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC133_Unique
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC133_intersect*
rm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC133_Union*

kmc_tools transform /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC133_Unique dump /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/CC133_MarkerKmers.txt
