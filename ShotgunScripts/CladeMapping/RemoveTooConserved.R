# Amy Campbell
# March 2023
# Remove any genes with >30X mean depth from S. epi or S. pettenkefori reads (2million 150bp simulated from each)


library(Biostrings)
library(stringr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Path to gene markers fasta
MarkerFastaPath=args[1]

# Path to the tsv containing covg stats by gene from S. epidermidis simulated reads
SEpiCovg=args[2]

# Path to the tsv containing covg stats by gene from S. pettenkefori simulated reads
SPetCovg=args[3]

RefGenome=args[4]


Markers = readDNAStringSet(MarkerFastaPath)
Sepiframe = read.table(SEpiCovg, col.names=c("identifer", "length", "mapped","placed", "min_cov", "max_cov", "mean_cov"))
SPetframe = read.table(SPetCovg, col.names=c("identifer", "length", "mapped","placed", "min_cov", "max_cov", "mean_cov"))

RemoveEpi = (Sepiframe %>% filter(mean_cov>30))$identifer
RemovePet = (SPetframe %>% filter(mean_cov>30))$identifer

Removes = c(RemoveEpi, RemovePet)
RemovesNames = paste0(Removes, " ", RefGenome)
KeepsNames = setdiff(names(Markers), RemovesNames)

FinalFastaPath=str_replace(MarkerFastaPath, ".fasta", "_Final.fasta")
NewMarkers=Markers[KeepsNames]

writeXStringSet(NewMarkers, file=FinalFastaPath)

