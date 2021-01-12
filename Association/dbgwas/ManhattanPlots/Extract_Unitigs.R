
# Amy Campbell
# 01/2021
# Extract the Unitigs and their associated IDs from dbGWAS/bugWAS into .fasta file
# so that we can map them back to a reference genome using BWA mem 

# Environment and required packages
###################################
setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")
library(dplyr)
library(tidyverse)
library(seqinr)
NodeInfoPath = "/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/Results/DBGWAS_Output_WithAnnotations/textualOutput/all_comps_nodes_info.tsv"

#################
NodeInfo = data.frame(read_tsv(NodeInfoPath))
NodeInfo = NodeInfo %>% filter(SequenceLength >20 )
OnlyEssential = NodeInfo %>% select(c(NodeId, Sequence))
write.fasta(as.list(OnlyEssential$Sequence), as.list(OnlyEssential$NodeId), file.out="data/Unitigs.fasta")

