# Amy Campbell
# December 2020
# Make metadata for input into DBGWAS 
bigframe = read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/BigFrame.csv")

setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")
phenotypes = read.csv("data/Phenotypes_01.15.20.csv")
tack_on = "/home/acampbe/FinalContigs/"
DBGWAS_table = phenotypes[c("DORN","XanthinPhenotype")]
DBGWAS_table = DBGWAS_table %>% filter(DORN != 2181)
DBGWAS_table = DBGWAS_table %>% filter(DORN != 685)
DBGWAS_table = DBGWAS_table %>% filter(DORN != 946)


DBGWAS_table$ID = paste0("DORN",DBGWAS_table$DORN)
DBGWAS_table$Phenotype = DBGWAS_table$XanthinPhenotype
DBGWAS_table$Path = paste0(tack_on, DBGWAS_table$ID)
DBGWAS_table$Path = paste0(DBGWAS_table$Path, "_cleaned.fasta")
DBGWAS_table = DBGWAS_table[c("ID", "Phenotype", "Path")]
write.table(DBGWAS_table, file="XanthinMapping", sep="\t",quote = F, row.names=F)

DBGWAS_table
View(DBGWAS_table)
DBGWAS_table
View(DBGWAS_table)

setdiff(bigframe$DORN, DBGWAS_table$ID)


