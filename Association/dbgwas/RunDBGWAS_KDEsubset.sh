# Amy Campbell
# June 2022
# Association testing with DBGWAS 
# Subsets based on KDE 

source /home/acampbe/software/miniconda3/bin/activate DBGWASEnv

# Staphyloxanthin
# 110 KDE-subset  isolates
# Zero-one-normalized
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains /home/acampbe/DFUStrainsWGS/mappings/XanthinZeroOneSubsetKDE.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Xanthin_kdeSubset/

# Siderophore
# Zero-one-normalized
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains /home/acampbe/DFUStrainsWGS/mappings/SiderophoreZeroOneSubsetKDE.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Siderophore_kdeSubset/

# Staphylokinase
# Just the 131 isolates containing sak (86 isolate subset of that)
# Zero-one normalized
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains /home/acampbe/DFUStrainsWGS/mappings/KinaseZeroOneSubsetKDE.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Kinase_kdeSubset/

# Biofilm
# Subset of 95 isolates(KDE-based)
# Zero-one-normalized
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains  /home/acampbe/DFUStrainsWGS/mappings/BiofilmZeroOneSubsetKDE.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Biofilm_Subset/

