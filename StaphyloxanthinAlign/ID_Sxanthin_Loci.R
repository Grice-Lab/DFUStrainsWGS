# Amy Campbell
# October 2020
# Seeing if we can pull the promoter out of the
# genome with everything else for the crtOPQMN operon

library(dplyr)

# setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/StaphyloxanthinAlign")

# Read in thse loci of the promoters
Promoters = read.table("AllPromoterSequences.tab", header=F)
colnames(Promoters) = c("Genome", "Gene", "Contig_Promoter", "PctID", "Length", "Start", "End", "queryCoverage","evalue", "bitscore")

# Read in the info about the rest of the genes in the cassette
XanthinGenes = read.csv("XanthinContigInfo.csv")
colnames(XanthinGenes) = c("Genome", "Samecontig", "Contig_Genes", "StartPositionGenes", "EndPositionGenes", "LengthGenes")

Merged = XanthinGenes %>% left_join(Promoters, by = "Genome")
Merged = Merged %>% mutate(AllSameContig = if_else(Contig_Promoter==Contig_Genes, TRUE, FALSE))

Merged = Merged %>% mutate(HighestLocus = max(c(Start,End, StartPositionGenes, EndPositionGenes)))
Merged_Select = Merged %>% select(c(Start,End, StartPositionGenes, EndPositionGenes))
Merged$HighestCoordinates = apply(Merged_Select, 1, max)
Merged$LowestCoordinates = apply(Merged_Select, 1, min)
Merged$TotalLengthChunk= Merged$HighestCoordinates - Merged$LowestCoordinates

# Add 20 bp on each end just in case :) 
Merged$StopOperon = Merged$HighestCoordinates 
Merged$StartOperon = Merged$LowestCoordinates
Merged$StopOperon20 = Merged$HighestCoordinates + 20
Merged$StartOperon20 = Merged$LowestCoordinates - 20


Merged_Output = Merged %>% select(c(Genome, Contig_Promoter, StartOperon, StopOperon ,StartOperon20,  StopOperon20))
write.csv(Merged_Output, file="ExtractionTable.csv")
