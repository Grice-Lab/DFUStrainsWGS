# Amy Campbell
# June 2022
# Association testing with DBGWAS 

source /home/acampbe/software/miniconda3/bin/activate DBGWASEnv

# Staphyloxanthin
# All 219 isolates
# Zero-one-normalized
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains  /home/acampbe/DFUStrainsWGS/mappings/TestZeroOneStaphyloxanthin.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Staphyloxanthin_Output_ZeroOne/

# Siderophore
# All 219 isolates
# Zero-one-normalized
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains  /home/acampbe/DFUStrainsWGS/mappings/TestZeroOneSiderophore.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Siderophore_Output_ZeroOne/

# Staphylokinase
# Just the 131 isolates containing sak (no further subsetting)
# Zero-one normalized (across the 219 isolates)
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains  /home/acampbe/DFUStrainsWGS/mappings/TestZeroOneStaphylokinase.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Kinase_Output_ZeroOne/

# Biofilm
# All 219
# Zero-one-normalized
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains  /home/acampbe/DFUStrainsWGS/mappings/TestZeroOneBiofilm.txt -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Biofilm_Output_ZeroOne/

