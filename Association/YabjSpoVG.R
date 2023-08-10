# Amy Campbell
# Looking at Yabj/SpoVG variation across the dataset

library(dplyr)
library(Biostrings)
library(stringr)
OperonSequenceJE2 = Biostrings::readDNAStringSet("/Users/amycampbell/Documents/DataInputGithub/data/RNASeqNE404/USA300_FPR3757.fasta", format="fasta")
Yabj_And_SpoVG = (OperonSequenceJE2$contig1)[532549:533303]

Yabj_And_SpoVG = Biostrings::DNAStringSet(Yabj_And_SpoVG)

Biostrings::writeXStringSet(c(Yabj_And_SpoVG),format="fasta" ,"/Users/amycampbell/Documents/DataInputGithub/data/RNASeqNE404/yabJ_spoVGregion_USA300_FPR3757.fasta")

Patient_Genome_Info = read.csv("~/Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")
Patient_Genome_Info$DORN = paste0("DORN", Patient_Genome_Info$Doern.lab.bank.)
Patient_Genome_Info = Patient_Genome_Info %>% select(patient_id, DORN) %>% unique()

roaryinfo= read.csv("Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv")

roaryinfoOperon = roaryinfo %>% filter(Gene %in% c("group_2993","extdb:pgaptmp_001576","spoVG","group_4414"))
roaryinfoOperonPresAbs = roaryinfoOperon
roaryinfoOperonPresAbs = roaryinfoOperonPresAbs %>% select(Gene,colnames(roaryinfoOperonPresAbs)[grepl(colnames(roaryinfoOperonPresAbs), pattern="DORN")])


JustPresAbs =roaryinfoOperonPresAbs[2:ncol(roaryinfoOperonPresAbs)]
JustPresAbs[JustPresAbs==""] <- 0
JustPresAbs[JustPresAbs!=0]  <- 1

JustPresAbs = JustPresAbs %>% mutate_all(function(x) as.numeric(as.character(x)))

GenomesMissingYabjSpoVG = colnames(JustPresAbs)[colSums(JustPresAbs)==0]

roaryinfoOperonPresAbs %>% filter(GenomesMissingYabjSpoVG)




JustGenomesWithoutYabjSpoVG = roaryinfo %>% select(GenomesMissingYabjSpoVG)


# Just spoVG
############
RoaryInfoSpoVG= roaryinfo %>% filter(Gene=="spoVG")
RoaryInfoSpoVG = RoaryInfoSpoVG %>% select(colnames(RoaryInfoSpoVG)[grepl("DORN",colnames(RoaryInfoSpoVG))])

SpoVGDf = data.frame(t(RoaryInfoSpoVG))
SpoVGDf$Genome = row.names(SpoVGDf)
colnames(SpoVGDf) = c("Mapping", "Genome")

# check if any of them have multiple mappings (none do)
SpoVGDf$containsTab = if_else(grepl(SpoVGDf$Mapping, pattern='\t'), "Yes", "No")
SpoVGDf$containsTab = NULL

SpoVGDf = SpoVGDf %>% filter(Mapping!="")
write.csv(SpoVGDf, row.names = F, quote = F, file="/Users/amycampbell/Documents/DataInputGithub/data/RNASeqNE404/SpoVGmappings.csv")

# Just yabJ
###########
RoaryInfoYabJ = roaryinfo %>% filter(Gene=="extdb:pgaptmp_001576")

RoaryInfoYabJ = RoaryInfoYabJ %>% select(colnames(RoaryInfoYabJ)[grepl("DORN",colnames(RoaryInfoYabJ))])

YabJDf = data.frame(t(RoaryInfoYabJ))
YabJDf$Genome = row.names(YabJDf)
colnames(YabJDf) = c("Mapping", "Genome")

# check if any of them have multiple mappings (none do)
YabJDf$containsTab = if_else(grepl(YabJDf$Mapping, pattern='\t'), "Yes", "No")
YabJDf$containsTab = NULL
YabJDf = YabJDf  %>% filter(Mapping!="")

write.csv(YabJDf, row.names = F, quote = F, file="/Users/amycampbell/Documents/DataInputGithub/data/RNASeqNE404/YabJmappings.csv")
