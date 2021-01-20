# Amy Campbell
# January 2021
# Extract unitig sequences from DBGWAS for alignment to 502A

source /home/acampbe/software/miniconda3/bin/activate DBGWASEnv


INPUT="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_31/textualOutput/all_comps_nodes_info.tsv"
OUTPUT="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_31/Unitigs.fasta"

Rscript --vanilla Extract_Unitigs.R --unitigs $INPUT --output $OUTPUT 
