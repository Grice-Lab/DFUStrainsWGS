# Amy Campbell
# January 2021
# Extract unitig sequences from DBGWAS for alignment to 502A

source /home/acampbe/software/miniconda3/bin/activate DBGWASEnv


INPUT_NODES="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_207_3-17/step1/graph.nodes"

INPUT_PATTERNS="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_207_3-17/step2/patterns.txt"

OUTPUT="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_207_3-17/Unitigs207_AllPatterns.fasta"


Rscript --vanilla Extract_Unitigs_AllPatterns.R --nodes $INPUT_NODES --patterns $INPUT_PATTERNS --output $OUTPUT 
