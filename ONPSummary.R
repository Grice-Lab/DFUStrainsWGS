# Amy Campbell
# July 2021
# Changing file names to DORNs and taking note of who got sequenced
# First step in assembling the hybrid assemblies
# Make mapping file to rename the files so that they correspond to DORN filenames 



library("dplyr")
library("stringr")
library("tidyr")
setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/")


All_ARMS = read.csv("CompletedLongReadSampleNames_CSV.csv", header=F)


All_ARMS = All_ARMS %>% filter(grepl(V1, pattern="ARM"))

ARMs_Sequenced = read.csv("ONPJulySequenceFiles.csv")
ARMs_Unsequenced = read.csv("ExcludedONPJulySequenceFiles.csv")

ARM_Mapping = read.csv("data/AmeliaStrains.csv")
ARM_Mapping$CorrespondingARM  = ARM_Mapping$Strain.Number

ARM_Mapping$CorrespondingARM = str_replace(ARM_Mapping$CorrespondingARM, "ARM0", "ARM")
ARM_Mapping = ARM_Mapping %>% select(CorrespondingARM, patient, doern.lab.bank..)

ARMs_Sequenced = ARMs_Sequenced %>% left_join(ARM_Mapping, by="CorrespondingARM")
ARMs_Sequenced = ARMs_Sequenced %>% drop_na()

ARMs_Sequenced$NewFileName = paste0("DORN", ARMs_Sequenced$doern.lab.bank..)
ARMs_Sequenced$DORN = ARMs_Sequenced$NewFileName
ARMs_Sequenced$NewFileName = paste0(ARMs_Sequenced$NewFileName, "_longread.fastq.gz")

ARMs_Unsequenced = ARMs_Unsequenced %>% left_join(ARM_Mapping, by="CorrespondingARM")

length(setdiff(unique(ARMs_Unsequenced$patient), unique(ARMs_Sequenced$patient)))

ARMs_Sequenced_FNames = ARMs_Sequenced %>% select(DORN, FileNameMIGS, NewFileName)
write.table(ARMs_Sequenced_FNames, row.names=F, col.names=F, quote=F, sep="\t", file="data/ONPFileMap.txt")

# 43 patients represented total
# 16 patients excluded total 
