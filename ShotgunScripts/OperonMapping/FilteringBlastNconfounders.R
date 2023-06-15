# Amy Campbell
# Have these long multifastas of sequences from NCBI that cover at least 150bp of the USA300 version of each sigB operon gene, 
# and have at least 70% identity to the sigB nt sequence 
# But want to filter them out based on their coming from species that are estimated to be in the metagenomes at all 
# 

library(stringr)
library(dplyr)
Metaphlan =read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Species_By_Sample_TimeSubj_MetaPhlan.csv")

Metaphlan$SubjectID=NULL
Metaphlan$Timepoint=NULL
Metaphlan$healed_by_12_weeks = NULL
Metaphlan$X=NULL

max_by_species = apply(Metaphlan,2,max)
max_by_species = max_by_species[max_by_species>.01]

SpeciesInclude = names(max_by_species)

# Filter RsbU homologs
######################
RsbU = Biostrings::readDNAStringSet("Documents/DataInputGithub/data/VariantsSigB/BlastNresults/rsbU_homologs_Blastn.fasta")

RsbUHeaders = names(RsbU)
RsbUNames = data.frame(Header=RsbUHeaders)

RsbUNames$Genus = sapply(RsbUNames$Header, function(x) str_split(x, " ")[[1]][2])
RsbUNames$Species= sapply(RsbUNames$Header, function(x) str_split(x, " ")[[1]][3])
RsbUNames$Metaphlan_Style = paste(RsbUNames$Genus, RsbUNames$Species, sep="_")
KeepRsbU = intersect(RsbUNames$Metaphlan_Style,SpeciesInclude )
RsbUNames  = RsbUNames %>% filter(Metaphlan_Style %in% KeepRsbU)

FilteredRsbUObj = RsbU[RsbUNames$Header]
Biostrings::writeXStringSet(FilteredRsbUObj, filepath="Documents/DataInputGithub/data/VariantsSigB/BlastNresults/FilteredRsbUHomologs.fasta", format = "fasta")

# Filter SigB homologs
######################
SigB = Biostrings::readDNAStringSet("Documents/DataInputGithub/data/VariantsSigB/BlastNresults/sigB_homologs_Blastn.fasta")
SigBHeaders = names(SigB)
SigBNames = data.frame(Header=SigBHeaders)

SigBNames$Genus = sapply(SigBNames$Header, function(x) str_split(x, " ")[[1]][2])
SigBNames$Species= sapply(SigBNames$Header, function(x) str_split(x, " ")[[1]][3])
SigBNames$Metaphlan_Style = paste(SigBNames$Genus, SigBNames$Species, sep="_")
KeepSigB = intersect(SigBNames$Metaphlan_Style,SpeciesInclude )
SigBNames  = SigBNames %>% filter(Metaphlan_Style %in% KeepSigB)
FilteredSigBObj = SigB[SigBNames$Header]
Biostrings::writeXStringSet(FilteredSigBObj, filepath="Documents/DataInputGithub/data/VariantsSigB/BlastNresults/FilteredSigBHomologs.fasta", format = "fasta")


# Filter RsbV homologs
######################
RsbV = Biostrings::readDNAStringSet("~/Documents/DataInputGithub/data/VariantsSigB/BlastNresults/rsbV_homologs_Blastn.fasta")
RsbVHeaders = names(RsbV)
RsbVNames = data.frame(Header=RsbVHeaders)

RsbVNames$Genus = sapply(RsbVNames$Header, function(x) str_split(x, " ")[[1]][2])
RsbVNames$Species= sapply(RsbVNames$Header, function(x) str_split(x, " ")[[1]][3])
RsbVNames$Metaphlan_Style = paste(RsbVNames$Genus, RsbVNames$Species, sep="_")
KeepRsbV = intersect(RsbVNames$Metaphlan_Style,SpeciesInclude )
RsbVNames = RsbVNames %>% filter(Metaphlan_Style %in% KeepRsbV)
FilteredRsbVObj = RsbV[RsbVNames$Header]
Biostrings::writeXStringSet(FilteredRsbVObj, filepath="Documents/DataInputGithub/data/VariantsSigB/BlastNresults/FilteredRsbVHomologs.fasta", format = "fasta")


# Filter RsbW homologs
######################
RsbW = Biostrings::readDNAStringSet("~/Documents/DataInputGithub/data/VariantsSigB/BlastNresults/rsbW_homologs_Blastn.fasta")
RsbWHeaders = names(RsbW)
RsbWNames = data.frame(Header=RsbWHeaders)

RsbWNames$Genus = sapply(RsbWNames$Header, function(x) str_split(x, " ")[[1]][2])
RsbWNames$Species= sapply(RsbWNames$Header, function(x) str_split(x, " ")[[1]][3])
RsbWNames$Metaphlan_Style = paste(RsbWNames$Genus, RsbWNames$Species, sep="_")
KeepRsbW = intersect(RsbWNames$Metaphlan_Style,SpeciesInclude )
RsbWNames = RsbWNames %>% filter(Metaphlan_Style %in% KeepRsbW)

FilteredRsbWObj = RsbW[RsbWNames$Header]
Biostrings::writeXStringSet(FilteredRsbWObj, filepath="Documents/DataInputGithub/data/VariantsSigB/BlastNresults/FilteredRsbWHomologs.fasta", format = "fasta")



