# Amy Campbell
# December 2020
# Make metadata for input into DBGWAS 
library(dplyr)
library(ggplot2)
setwd("~/Desktop/GriceLabGit/DFUStrainsWGS")
bigframe = read.csv("data/BigFrame.csv")

setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")
phenotypes = read.csv("data/Phenotypes_01.15.20.csv")
tack_on = "/home/acampbe/FinalContigs/"

DBGWAS_table = phenotypes[c("DORN","XanthinPhenotype")]
DBGWAS_table = DBGWAS_table %>% filter(DORN != 2181)

# Not S. aureus
DBGWAS_table = DBGWAS_table %>% filter(DORN != 685)
DBGWAS_table = DBGWAS_table %>% filter(DORN != 946)

DBGWAS_table$ID = paste0("DORN",DBGWAS_table$DORN)
DBGWAS_table$Phenotype = DBGWAS_table$XanthinPhenotype
DBGWAS_table$Path = paste0(tack_on, DBGWAS_table$ID)
DBGWAS_table$Path = paste0(DBGWAS_table$Path, "_cleaned.fasta")
DBGWAS_table = DBGWAS_table[c("ID", "Phenotype", "Path")]
write.table(DBGWAS_table, file="XanthinMapping", sep="\t",quote = F, row.names=F)

# Amelia working on finishing these two 
setdiff(bigframe$DORN, DBGWAS_table$ID)

ggplot(DBGWAS_table, aes(x=Phenotype)) + geom_histogram(fill="darkgoldenrod", color="black") + xlab("Staphyloxanthin Production (Absorbance, Normalized to SA502A)") + ggtitle("Staphyloxanthin Production in 221 S. aureus Isolates")

phenotypes_Updated = read.csv("data/staphyloxanthin_averages_1.13.21.csv")
colnames(phenotypes_Updated) = c("ARM", "DORN", "SUBJ","TIME", "AVERAGE", "SD")
DBGWAS_table_updated = phenotypes_Updated[c("DORN","AVERAGE")]
DBGWAS_table_updated = DBGWAS_table_updated %>% filter(DORN != 2181)

# Not S. aureus
DBGWAS_table_updated = DBGWAS_table_updated %>% filter(DORN != 685)
DBGWAS_table_updated = DBGWAS_table_updated %>% filter(DORN != 946)
ggplot(DBGWAS_table_updated, aes(x=Phenotype)) + geom_histogram(fill="darkgoldenrod", color="black") + xlab("Staphyloxanthin Production (Absorbance, Normalized to SA502A)") + ggtitle("Staphyloxanthin Production in 219 S. aureus Isolates")


DBGWAS_table_updated$ID = paste0("DORN",DBGWAS_table_updated$DORN)
DBGWAS_table_updated$Phenotype = DBGWAS_table_updated$AVERAGE
DBGWAS_table_updated$Path = paste0(tack_on, DBGWAS_table_updated$ID)
DBGWAS_table_updated$Path = paste0(DBGWAS_table_updated$Path, "_cleaned.fasta")
DBGWAS_table_updated = DBGWAS_table_updated[c("ID", "Phenotype", "Path")]
write.table(DBGWAS_table_updated, file="XanthinMappingFull", sep="\t",quote = F, row.names=F)

# Update 02/2021 With removed 14 isolates (contaminated)
########################################################
# Also, log-transforming the xanthin phenotype for normality. not SURE if this matters? Should try first running it on what I was using before. 

Contaminated = read.csv("data/ContaminatedIsolates2-12-21.csv")
DBGWAS_table_Remove14 = DBGWAS_table_updated %>% filter(!(ID %in% Contaminated$DORN))
ggplot(DBGWAS_table_Remove14, aes(x=as.numeric(as.character(Phenotype)))) + geom_histogram(fill="darkgoldenrod", color="black") + xlab("Staphyloxanthin Production (Absorbance, Normalized to SA502A)") + ggtitle("Staphyloxanthin Production in 219 S. aureus Isolates")
ggplot(DBGWAS_table_Remove14, aes(x=log(as.numeric(as.character(Phenotype))))) + geom_histogram(fill="darkgoldenrod", color="black") + xlab("Log-transformed Staphyloxanthin Production (Absorbance, Normalized to SA502A)") + ggtitle("Staphyloxanthin Production in 219 S. aureus Isolates")

DBGWAS_table_Remove14_logtransformed = DBGWAS_table_Remove14 
DBGWAS_table_Remove14_logtransformed$Phenotype = as.numeric(as.character(DBGWAS_table_Remove14_logtransformed$Phenotype))
DBGWAS_table_Remove14_logtransformed$Phenotype = log(DBGWAS_table_Remove14_logtransformed$Phenotype)
write.table(DBGWAS_table_Remove14_logtransformed, file="XanthinMapping207Log", sep="\t",quote = F, row.names=F)
write.table(DBGWAS_table_Remove14,file="XanthinMapping207", sep="\t",quote = F, row.names=F )

DBGWAS_table_Remove14$Phenotype = as.numeric(as.character(DBGWAS_table_Remove14$Phenotype))

