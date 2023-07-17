# Amy Campbell
# Having summarized the presence/absence of genes by patient, are there any genes whose presence in any isolate of a patient is associated with healing or nonhealing?
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

# Filter to genes in  more than 10% but less than 90%
GeneByPatientDF = GeneByPatientDF %>% filter(RowSum <= 54 & RowSum >= 6)

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
