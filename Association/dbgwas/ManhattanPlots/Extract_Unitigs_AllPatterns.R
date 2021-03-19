
# Amy Campbell
# 01/2021
# Extract the Unitigs and their associated IDs from dbGWAS/bugWAS into .fasta file
# so that we can map them back to a reference genome using BWA mem 

# Environment and required packages
###################################

library(dplyr)
library(tidyverse)
library(seqinr)
library(optparse)


option_list = list(make_option(c("-n", "--nodes"), type="character", help="Absolute path to graph.nodes file from DBGWAS"), make_option(c("-p", "--patterns"), type="character", help="Absolute path to patterns.txt file from DBGWAS bugWAS step"), make_option(c("-o", "--output"), type="character", help="Absolute path to desired output file"))


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

NodeInfoPath = opt$nodes

OutputPath= opt$output

PatternsPath= opt$patterns

#################

bugWASpatterns =  read_delim(PatternsPath, " ", col_names=T)
bugWASpatterns = bugWASpatterns$pattern

nodes = read_tsv(NodeInfoPath, col_names=F)
colnames(nodes) = c("pattern", "sequence")
nodes$sequencelength = nchar(nodes$sequence)
nodes = nodes %>% filter(sequencelength>20)
nodes = nodes %>% filter (pattern %in% bugWASpatterns)
nodes$NodeId = paste0("n", nodes$pattern)

OnlyEssential = nodes %>% select(c(NodeId, sequence))


write.fasta(as.list(OnlyEssential$sequence), as.list(OnlyEssential$NodeId), file.out=OutputPath)

