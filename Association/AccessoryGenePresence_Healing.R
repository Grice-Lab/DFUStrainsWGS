# Amy Campbell
# Having summarized the presence/absence of genes by patient, are there any genes whose presence in any isolate of a patient is associated with healing or nonhealing?
library(stringr)
library(dplyr)

set.seed(19104)
UpdatedPhenotypes=read.csv("~/Desktop/GriceLabGit/Staphyloxanthin/Data/InVitroData/Phenotypes_Data.csv")

UpdatedPhenotypes = UpdatedPhenotypes %>% select(patient, week_healed) %>% unique()

UpdatedPhenotypes$HealedBy12= if_else(!is.na(UpdatedPhenotypes$week_healed) & UpdatedPhenotypes$week_healed <=12, "Yes", "No")
PatientHealMap = UpdatedPhenotypes %>% select(patient, HealedBy12)
PatientHealMap$patient = sapply(PatientHealMap$patient, as.character)

GeneByPatientDF =read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient.csv")
cols_select = colnames(GeneByPatientDF)
cols_select = cols_select[grepl("patient", cols_select)]
GeneByPatientDF_Clean = GeneByPatientDF %>% select(cols_select)
saveGenes =GeneByPatientDF$Gene
savePatients=colnames(GeneByPatientDF_Clean)

GeneByPatientDF$RowSum = rowSums(GeneByPatientDF[5:ncol(GeneByPatientDF)])

# Filter to genes in  at least 10% but no more than 90%
GeneByPatientDF = GeneByPatientDF %>% filter(RowSum <= 60*.9 & RowSum >= 60*.1)

GeneByPatientDF$RowSum=NULL

JustGenes=as.matrix(GeneByPatientDF[,5:ncol(GeneByPatientDF)])
RowsToCheck = 1:nrow(JustGenes)

NewMatrix = c()

for(rownum in 1:nrow(JustGenes)){
  print(rownum)
  # If I haven't already collapsed identically distributed rows 
  if(rownum %in% RowsToCheck ){
    
    row_check=JustGenes[rownum,]
    NumIdentical=which(apply(JustGenes, 1, function(x) identical(x,row_check)))
    if(length(NumIdentical) >1){
      
      # Paste together gene names with "-" between them
      GroupGeneName=paste(GeneByPatientDF[NumIdentical,"Gene"], collapse="-")
      
      # Remove all the identical rows bc now theyre represented
      RowsToCheck = RowsToCheck[ !(RowsToCheck %in% NumIdentical)]
      
      # Add this to the matrix of gene-by-patient "haplotypes"
      NewMatrix = rbind(NewMatrix, c(GroupGeneName,row_check ))
    }else{
      RowsToCheck = RowsToCheck[RowsToCheck!=rownum]
      GeneName=GeneByPatientDF[rownum,"Gene"]
      NewMatrix = rbind(NewMatrix, c(GeneName,row_check ))
    }

  }
 
}


NewMatrix=data.frame(NewMatrix)
haplotypeNames=NewMatrix$V1

transposedDF = data.frame(t(NewMatrix[,2:ncol(NewMatrix)]))

colnames(transposedDF) = haplotypeNames

transposedDF$patient = sapply(row.names(transposedDF), function(x) str_split(x,"_")[[1]][2])
transposedDF =transposedDF %>% left_join(PatientHealMap,by="patient")


colnames(transposedDF)[ncol(transposedDF)-2]

genesTested = c()
chisqP=c()
cols_genes = colnames(transposedDF)
for(indexitem in 1:(ncol(transposedDF) -2 )){
  summary_dist=table(transposedDF[c(indexitem,ncol(transposedDF))])
  ChiSqResult = chisq.test(summary_dist, simulate.p.value = T)
  genesTested = append( genesTested, cols_genes[indexitem])
  chisqP = append(chisqP, ChiSqResult$p.value)
}
ChisqP_DF = data.frame(Gene=genesTested, pvalue=chisqP)

ChisqP_DF$Padjusted = p.adjust(ChisqP_DF$pvalue, length(ChisqP_DF$pvalue), method="BH")




# Now do JUST CC5-associated genes by patient
##############################################

GeneByPatientDF_CC5 =read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient_CC5.csv")
cols_selectCC5 = colnames(GeneByPatientDF_CC5)
cols_selectCC5 = cols_selectCC5[grepl("patient", cols_selectCC5)]
GeneByPatientDF_Clean_CC5 = GeneByPatientDF_CC5 %>% select(cols_selectCC5)
saveGenesCC5 =GeneByPatientDF_CC5$Gene
savePatientsCC5=colnames(GeneByPatientDF_Clean_CC5)


# "haplotype" filtering
GeneByPatientDF_CC5$RowSum = rowSums(GeneByPatientDF_CC5[5:ncol(GeneByPatientDF_CC5)])

# Filter to genes in  more than 10% but less than 90%
GeneByPatientDF_CC5 = GeneByPatientDF_CC5 %>% filter(RowSum <= (22*.8) & RowSum >= (22*.2))

GeneByPatientDF_CC5$RowSum=NULL

JustGenes=as.matrix(GeneByPatientDF_CC5[,5:ncol(GeneByPatientDF_CC5)])
RowsToCheck = 1:nrow(JustGenes)

NewMatrix = c()

for(rownum in 1:nrow(JustGenes)){
  print(rownum)
  # If I haven't already collapsed identically distributed rows 
  if(rownum %in% RowsToCheck ){
    
    row_check=JustGenes[rownum,]
    NumIdentical=which(apply(JustGenes, 1, function(x) identical(x,row_check)))
    if(length(NumIdentical) >1){
      
      # Paste together gene names with "-" between them
      GroupGeneName=paste(GeneByPatientDF_CC5[NumIdentical,"Gene"], collapse="-")
      
      # Remove all the identical rows bc now theyre represented
      RowsToCheck = RowsToCheck[ !(RowsToCheck %in% NumIdentical)]
      
      # Add this to the matrix of gene-by-patient "haplotypes"
      NewMatrix = rbind(NewMatrix, c(GroupGeneName,row_check ))
    }else{
      RowsToCheck = RowsToCheck[RowsToCheck!=rownum]
      GeneName=GeneByPatientDF_CC5[rownum,"Gene"]
      NewMatrix = rbind(NewMatrix, c(GeneName,row_check ))
    }
  }
}

NewMatrix=data.frame(NewMatrix)
haplotypeNames=NewMatrix$V1

transposedDF = data.frame(t(NewMatrix[,2:ncol(NewMatrix)]))

colnames(transposedDF) = haplotypeNames

transposedDF$patient = sapply(row.names(transposedDF), function(x) str_split(x,"_")[[1]][2])
transposedDF =transposedDF %>% left_join(PatientHealMap,by="patient")


colnames(transposedDF)[ncol(transposedDF)-2]

genesTested = c()
chisqP=c()
cols_genes = colnames(transposedDF)
for(indexitem in 1:(ncol(transposedDF) -2 )){
  summary_dist=table(transposedDF[c(indexitem,ncol(transposedDF))])
  ChiSqResult = chisq.test(summary_dist, simulate.p.value = T)
  genesTested = append( genesTested, cols_genes[indexitem])
  chisqP = append(chisqP, ChiSqResult$p.value)
}
ChisqP_DF = data.frame(Gene=genesTested, pvalue=chisqP)

ChisqP_DF$Padjusted = p.adjust(ChisqP_DF$pvalue, length(ChisqP_DF$pvalue), method="BH")


# Now just CC1
#############################################################

GeneByPatientDF_CC1 =read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient_CC1.csv")
cols_selectCC1 = colnames(GeneByPatientDF_CC1)
cols_selectCC1 = cols_selectCC1[grepl("patient", cols_selectCC1)]
GeneByPatientDF_Clean_CC1 = GeneByPatientDF_CC1 %>% select(cols_selectCC1)
saveGenesCC1 =GeneByPatientDF_CC1$Gene
savePatientsCC1=colnames(GeneByPatientDF_Clean_CC1)


# "haplotype" filtering
GeneByPatientDF_CC1$RowSum = rowSums(GeneByPatientDF_CC1[5:ncol(GeneByPatientDF_CC1)])

# Filter to genes in  more than 10% but less than 90%
GeneByPatientDF_CC1 = GeneByPatientDF_CC1 %>% filter(RowSum <= (11*.8) & RowSum >= (11*.2))

GeneByPatientDF_CC1$RowSum=NULL

JustGenes=as.matrix(GeneByPatientDF_CC1[,5:ncol(GeneByPatientDF_CC1)])
RowsToCheck = 1:nrow(JustGenes)

NewMatrix = c()

for(rownum in 1:nrow(JustGenes)){
  print(rownum)
  # If I haven't already collapsed identically distributed rows 
  if(rownum %in% RowsToCheck ){
    
    row_check=JustGenes[rownum,]
    NumIdentical=which(apply(JustGenes, 1, function(x) identical(x,row_check)))
    if(length(NumIdentical) >1){
      
      # Paste together gene names with "-" between them
      GroupGeneName=paste(GeneByPatientDF_CC1[NumIdentical,"Gene"], collapse="-")
      
      # Remove all the identical rows bc now theyre represented
      RowsToCheck = RowsToCheck[ !(RowsToCheck %in% NumIdentical)]
      
      # Add this to the matrix of gene-by-patient "haplotypes"
      NewMatrix = rbind(NewMatrix, c(GroupGeneName,row_check ))
    }else{
      RowsToCheck = RowsToCheck[RowsToCheck!=rownum]
      GeneName=GeneByPatientDF_CC1[rownum,"Gene"]
      NewMatrix = rbind(NewMatrix, c(GeneName,row_check ))
    }
  }
}

NewMatrix=data.frame(NewMatrix)
haplotypeNames=NewMatrix$V1

transposedDF = data.frame(t(NewMatrix[,2:ncol(NewMatrix)]))

colnames(transposedDF) = haplotypeNames

transposedDF$patient = sapply(row.names(transposedDF), function(x) str_split(x,"_")[[1]][2])
transposedDF =transposedDF %>% left_join(PatientHealMap,by="patient")


colnames(transposedDF)[ncol(transposedDF)-2]

genesTested = c()
chisqP=c()
cols_genes = colnames(transposedDF)
for(indexitem in 1:(ncol(transposedDF) -2 )){
  summary_dist=table(transposedDF[c(indexitem,ncol(transposedDF))])
  ChiSqResult = chisq.test(summary_dist, simulate.p.value = T)
  genesTested = append( genesTested, cols_genes[indexitem])
  chisqP = append(chisqP, ChiSqResult$p.value)
}
ChisqP_DF = data.frame(Gene=genesTested, pvalue=chisqP)

ChisqP_DF$Padjusted = p.adjust(ChisqP_DF$pvalue, length(ChisqP_DF$pvalue), method="BH")



# Now just CC8
#############################################################

GeneByPatientDF_CC8 =read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient_CC8.csv")
cols_selectCC8 = colnames(GeneByPatientDF_CC8)
cols_selectCC8 = cols_selectCC8[grepl("patient", cols_selectCC8)]
GeneByPatientDF_Clean_CC8 = GeneByPatientDF_CC8 %>% select(cols_selectCC8)
saveGenesCC8 =GeneByPatientDF_CC8$Gene
savePatientsCC8=colnames(GeneByPatientDF_Clean_CC8)


# "haplotype" filtering
GeneByPatientDF_CC8$RowSum = rowSums(GeneByPatientDF_CC8[5:ncol(GeneByPatientDF_CC8)])

# Filter to genes in  more than 10% but less than 90%
GeneByPatientDF_CC8 = GeneByPatientDF_CC8 %>% filter(RowSum <= (10*.8) & RowSum >= (10*.20))

GeneByPatientDF_CC8$RowSum=NULL

JustGenes=as.matrix(GeneByPatientDF_CC8[,5:ncol(GeneByPatientDF_CC8)])
RowsToCheck = 1:nrow(JustGenes)

NewMatrix = c()

for(rownum in 1:nrow(JustGenes)){
  print(rownum)
  # If I haven't already collapsed identically distributed rows 
  if(rownum %in% RowsToCheck ){
    
    row_check=JustGenes[rownum,]
    NumIdentical=which(apply(JustGenes, 1, function(x) identical(x,row_check)))
    if(length(NumIdentical) >1){
      
      # Paste together gene names with "-" between them
      GroupGeneName=paste(GeneByPatientDF_CC8[NumIdentical,"Gene"], collapse="-")
      
      # Remove all the identical rows bc now theyre represented
      RowsToCheck = RowsToCheck[ !(RowsToCheck %in% NumIdentical)]
      
      # Add this to the matrix of gene-by-patient "haplotypes"
      NewMatrix = rbind(NewMatrix, c(GroupGeneName,row_check ))
    }else{
      RowsToCheck = RowsToCheck[RowsToCheck!=rownum]
      GeneName=GeneByPatientDF_CC8[rownum,"Gene"]
      NewMatrix = rbind(NewMatrix, c(GeneName,row_check ))
    }
  }
}

NewMatrix=data.frame(NewMatrix)
haplotypeNames=NewMatrix$V1

transposedDF = data.frame(t(NewMatrix[,2:ncol(NewMatrix)]))

colnames(transposedDF) = haplotypeNames

transposedDF$patient = sapply(row.names(transposedDF), function(x) str_split(x,"_")[[1]][2])
transposedDF =transposedDF %>% left_join(PatientHealMap,by="patient")


colnames(transposedDF)[ncol(transposedDF)-2]

genesTested = c()
chisqP=c()
cols_genes = colnames(transposedDF)
for(indexitem in 1:(ncol(transposedDF) -2 )){
  summary_dist=table(transposedDF[c(indexitem,ncol(transposedDF))])
  ChiSqResult = chisq.test(summary_dist, simulate.p.value = T)
  genesTested = append( genesTested, cols_genes[indexitem])
  chisqP = append(chisqP, ChiSqResult$p.value)
}
ChisqP_DF = data.frame(Gene=genesTested, pvalue=chisqP)

ChisqP_DF$Padjusted = p.adjust(ChisqP_DF$pvalue, length(ChisqP_DF$pvalue), method="BH")



