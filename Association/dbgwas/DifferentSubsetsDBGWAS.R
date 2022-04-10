# Amy Campbell
# April 2022
# Writing out unitigs with the lowest 20 unique p values for each subset I tested dbgwas on 

Subset98 = read.csv2("Documents/DFUData/GWAS/all_comps_nodes_info98.tsv",sep='\t',row.names=as.character(1:7683))
Subset104 = read.csv2("Documents/DFUData/GWAS/all_comps_nodes_info104.tsv",sep='\t',row.names=as.character(1:10618))
Subset219 = read.csv2("Documents/DFUData/GWAS/all_comps_nodes_info219.tsv",sep='\t',row.names=as.character(1:3169))

ColNamesVect = c("NodeID", "AlleleFreq", "Pheno0Count","Pheno0TotalCount","Pheno1Count", "Pheno1TotalCount","NACount", "NATotalCount", "Significant", "pvalue", "qValue", "EstEffect", "WaldStat", "Sequence", "SequenceLength", "Annotations")


Subset98 = Subset98[,2:17]
Subset104 = Subset104[,2:17]
Subset219 = Subset219[,2:17]


colnames(Subset98)= ColNamesVect
colnames(Subset104)= ColNamesVect
colnames(Subset219)= ColNamesVect


Subset98$pvalue = sapply(Subset98$pvalue, function(x) as.numeric(as.character(x)))
Subset104$pvalue = sapply(Subset104$pvalue, function(x) as.numeric(as.character(x)))
Subset219$pvalue = sapply(Subset219$pvalue, function(x) as.numeric(as.character(x)))

Subset98$qValue = sapply(Subset98$qValue, function(x) as.numeric(as.character(x)))
Subset104$qValue = sapply(Subset104$qValue, function(x) as.numeric(as.character(x)))
Subset219$qValue = sapply(Subset219$qValue, function(x) as.numeric(as.character(x)))


Subset98LowestPs = sort(unique(Subset98$pvalue))[1:20]
Subset104LowestPs = sort(unique(Subset104$pvalue))[1:20]
Subset219LowestPs =  sort(unique(Subset219$pvalue))[1:20]

Subset98Smallest = Subset98 %>% filter(pvalue <= Subset98LowestPs[20])
Subset104Smallest = Subset104 %>% filter(pvalue <= Subset104LowestPs[20])
Subset219Smallest =  Subset219 %>% filter(pvalue <= Subset219LowestPs[20])


WriteOutTopUnitigs98 = Subset98Smallest %>% select(NodeID, Sequence)

WriteOutTopUnitigs104 = Subset104Smallest %>% select(NodeID, Sequence)

WriteOutTopUnitigs219 = Subset219Smallest %>% select(NodeID, Sequence)

Open98= file("TopUnitigs98.fasta")

for(r in 1:nrow(WriteOutTopUnitigs98)){
  print(WriteOutTopUnitigs98[r, "Sequence"])
  write(c(paste0(">",WriteOutTopUnitigs98[r, "NodeID"] ), WriteOutTopUnitigs98[r, "Sequence"]), file="TopUnitigs98.fasta" , append=T)
}


Open104= file("TopUnitigs104.fasta")
for(r in 1:nrow(WriteOutTopUnitigs104)){
  print(WriteOutTopUnitigs104[r, "Sequence"])
  write(c(paste0(">",WriteOutTopUnitigs104[r, "NodeID"] ), WriteOutTopUnitigs104[r, "Sequence"]), file="TopUnitigs104.fasta" , append=T)
}


Open219= file("TopUnitigs219.fasta")
for(r in 1:nrow(WriteOutTopUnitigs219)){
  print(WriteOutTopUnitigs219[r, "Sequence"])
  write(c(paste0(">",WriteOutTopUnitigs219[r, "NodeID"] ), WriteOutTopUnitigs219[r, "Sequence"]), file="TopUnitigs219.fasta" , append=T)
}

Shared =intersect(Subset104Smallest$Sequence,Subset219Smallest$Sequence )
intersect(Shared, Subset98Smallest$Sequence)
intersect(Subset219Smallest$Sequence , Subset104Smallest$Sequence)



