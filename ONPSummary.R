# Amy Campbell
# July 2021
# Changing file names to DORNs and taking note of who got sequenced
# First step in assembling the hybrid assemblies
# Make mapping file to rename the files so that they correspond to DORN filenames 



library("dplyr")
library("stringr")
library("tidyr")
library("ggplot2")
setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/")


All_ARMS = read.csv("CompletedLongReadSampleNames_CSV.csv", header=F)

ShouldHaveExcluded = read.csv("data/ContaminatedIsolates2-12-21.csv")

All_ARMS = All_ARMS %>% filter(grepl(V1, pattern="ARM"))
colnames(All_ARMS) = c("CorrespondingARM","Sequenced")

ARMs_Sequenced = read.csv("ONPJulySequenceFiles.csv")
ARMs_Unsequenced = read.csv("ExcludedONPJulySequenceFiles.csv")


ARM_Mapping = read.csv("data/AmeliaStrains.csv")
ARM_Mapping$CorrespondingARM  = ARM_Mapping$Strain.Number

ARM_Mapping$CorrespondingARM = str_replace(ARM_Mapping$CorrespondingARM, "ARM0", "ARM")
ARM_Mapping = ARM_Mapping %>% select(CorrespondingARM, patient, doern.lab.bank..)

All_ARMS = All_ARMS %>% left_join(ARM_Mapping, by="CorrespondingARM")
All_ARMS %>% filter( doern.lab.bank.. %in% ShouldHaveExcluded$Doern.lab.bank.)



ARMs_Sequenced = ARMs_Sequenced %>% left_join(ARM_Mapping, by="CorrespondingARM")
ARMs_Sequenced = ARMs_Sequenced %>% drop_na()

ARMs_Sequenced$NewFileName = paste0("DOR rep(tab$value, tab$freq)N", ARMs_Sequenced$doern.lab.bank..)
ARMs_Sequenced$DORN = ARMs_Sequenced$NewFileName
ARMs_Sequenced$NewFileName = paste0(ARMs_Sequenced$NewFileName, "_longread.fastq")
ARMs_Sequenced$FileNameMIGS = str_replace(ARMs_Sequenced$FileNameMIGS, "fastq.gz", "fastq")
ARMs_Unsequenced = ARMs_Unsequenced %>% left_join(ARM_Mapping, by="CorrespondingARM")

length(setdiff(unique(ARMs_Unsequenced$patient), unique(ARMs_Sequenced$patient)))

ARMs_Sequenced_FNames = ARMs_Sequenced %>% select(DORN, FileNameMIGS, NewFileName)
write.table(ARMs_Sequenced_FNames, row.names=F, col.names=F, quote=F, sep="\t", file="data/ONPFileMap.txt")



write.csv(ARMs_Unsequenced, file="MigsUnsequencedARMS.csv", quote=F, row.names=F, col.names=T)


ReadCounts = read.table("data/AllSampleReads.tsv")
colnames(ReadCounts) = c("Count", "Length", "DORN")
ReadCounts$Index = row.names(ReadCounts)
ggplot(ReadCounts, aes(x=Index, y=))

ReadCountsDF = data.frame(rep(ReadCounts$Length, ReadCounts$Count))
colnames(ReadCountsDF) = c("Length", "DORN")

p <- ggplot(ReadCountsDF, aes(x=DORN, y=log10(Length))) + geom_boxplot(fill="khaki1") + theme_classic() + scale_y_continuous(breaks=0:7) + ggtitle("Long Read Lengths") + theme(axis.text.x=element_text(angle=90)) + geom_line(y=3, linetype="dashed")
p
ggplot(ReadCountsDF, aes(x=Length))


ARMs_Sequenced %>% filter(patient==141)
# 43 patients represented total
# 16 patients excluded total 
All_ARMS %>% filter( doern.lab.bank.. %in% ShouldHaveExcluded$Doern.lab.bank.)
ARMs_Unsequenced %>% filter( doern.lab.bank.. %in% ShouldHaveExcluded$Doern.lab.bank.)



NewlySequencedSeptember = read.csv2("data/SecondONPrun.txt",header=F)
NewlySequencedSeptember$V1 = sapply(NewlySequencedSeptember, as.character)
NewlySequencedSeptember$CorrespondingARM = sapply(NewlySequencedSeptember$V1, function(x) (strsplit(x,"_"))[[1]][1])
NewlySequencedSeptember = NewlySequencedSeptember %>% left_join(ARM_Mapping, by="CorrespondingARM")
NewlySequencedSeptember$DORN = paste0("DORN", NewlySequencedSeptember$doern.lab.bank..)
NewlySequencedSeptember$NewFileName = paste0(NewlySequencedSeptember$DORN, "_longread.fastq.gz")

NewlySequencedSeptember = NewlySequencedSeptember %>% select(DORN, V1, NewFileName)
write.table(NewlySequencedSeptember, row.names=F, col.names=F, quote=F, sep="\t", file="data/ONPFileMapRun2.txt")

