# Amy Campbell
# December 2020
# Make metadata for input into DBGWAS 

setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")
phenotypes = read.csv("data/phenotype_variation_11.06.20.csv")

tack_on = "/home/acampbe/FinalContigs/"
DBGWAS_table = phenotypes[c("DORN","xanthin")]
DBGWAS_table$ID = paste0("DORN",DBGWAS_table$DORN)
DBGWAS_table$Phenotype = DBGWAS_table$xanthin
DBGWAS_table$Path = paste0(tack_on, DBGWAS_table$ID)
DBGWAS_table$Path = paste0(DBGWAS_table$Path, "_cleaned.fasta")
DBGWAS_table = DBGWAS_table[c("ID", "Phenotype", "Path")]
write.table(DBGWAS_table, file="XanthinMapping", sep="\t",quote = F, row.names=F)

