# Amy Campbell 2023
# Analyzing variants found in metagenomic shotgun reads in comparisons between early and late stage timepoints of DFU samples
# within the SigB operon 
library(dplyr)
library(stringr)

comparisonpath="/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/VariantComparisons/"

AllFiles = list.files(comparisonpath)

AdequateDepth = Filter(function(x) grepl("_Final.tsv",x), AllFiles)

GetComplement = function(string){
  print(string)
  if(grepl("C", string)){ 
  string = stringr::str_replace(string=string,pattern="C", replacement="cytosine")}
  if(grepl("G", string)){ 
    string = stringr::str_replace(string=string,pattern="G", replacement="guanine")}
  if(grepl("T", string)){ 
    string = stringr::str_replace(string=string,pattern="T", replacement="thymine")}
  if(grepl("A", string)){ 
    string = stringr::str_replace(string=string,pattern="A", replacement="adenine")}
  
  string = stringr::str_replace(string=string,pattern="cytosine",replacement="G")
  string = stringr::str_replace(string=string,pattern="guanine",replacement="C")
  string = stringr::str_replace(string=string,pattern="thymine",replacement="A")
  string = stringr::str_replace(string=string,pattern="adenine",replacement="T")
  print(string)
  return(string)
}

FilterConfounders = function(AltMetagenome, AltConfounders){
  Metagenomelist=str_split(AltMetagenome, pattern=",")[[1]]
  confounderslist=str_split(AltConfounders, pattern=",")[[1]]
  filteredList=setdiff(Metagenomelist,confounderslist )
  newalt = paste(filteredList, collapse=",")
  return(newalt)

}

Comparisons_sigB = Filter(function(x) grepl("sigB",x), AdequateDepth)
SigBcompareDF = data.frame(FileName=Comparisons_sigB)
SigBcompareDF$Patient = sapply(SigBcompareDF$FileName, function(x) str_split(x, "_")[[1]][3])
SigBcompareDF$Comparison = sapply(SigBcompareDF$FileName, function(x) paste((str_split(x, "_")[[1]][4]),(str_split(x, "_")[[1]][5]),  sep=":"))
table(SigBcompareDF$Comparison)

#Comparisons_rsbU = Filter(function(x) grepl("rsbU",x), AdequateDepth)
#rsbUcompareDF = data.frame(FileName=Comparisons_rsbU)


# SigB
#######

Confounders_SigB =read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/ConfoundersBlastN_Filtered/sigBsnps", sep="\t", comment.char = "#", header=F) #%>% select(V1,V2,V3,V4,V5,V6,V7,V8,V9)
Presence = Confounders_SigB[, 10:ncol(Confounders_SigB)]
Presence[Presence>0] <- 1
TotalPresent = rowSums(Presence)
Confounders_SigB = Confounders_SigB %>% select(V1,V2,V3,V4,V5,V6,V7,V8,V9)
Confounders_SigB$TotalPresent = TotalPresent

colnames(Confounders_SigB) = c("Chromosome", "Position", "ID","REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TotalPresent")

Confounders_SigB$BasePosition=2185766-(Confounders_SigB$Position - 1)
  
Comparisons_sigB = Filter(function(x) grepl("sigB",x), AdequateDepth)

SigBcompareDF$Patient = sapply(SigBcompareDF$FileName, function(x) str_split(x, "_")[[1]][3])
unique(SigBcompareDF$Patient)

SigBVariantsDF = data.frame()
for(patient in unique(SigBcompareDF$Patient)){
  subDF = SigBcompareDF %>% filter(Patient==patient)
  for(i in 1:nrow(subDF)){
    fname=(subDF$FileName)[i]
    ReadDF = read.table( paste0(comparisonpath, (subDF$FileName)[i]))
    comparison= paste((str_split(fname, "_")[[1]][4]),(str_split(fname, "_")[[1]][5]),  sep=":")   
    ReadDF$Patient = patient
    ReadDF$Comparison= comparison
    SigBVariantsDF = rbind(SigBVariantsDF,ReadDF )
    
  }
  
}

SigBVariantsDF$TotalAltDepth = sapply(SigBVariantsDF$I16, function(x) sum(sapply(str_split(x, ",")[[1]][c(3,4)], as.numeric) ))

SigBVariantsDFModerate_High = SigBVariantsDF %>% filter(PredictedEffect!="LOW" )

SigBVariantsDFHigh = SigBVariantsDF %>% filter(PredictedEffect=="HIGH" & TotalAltDepth>=2 )


SigBVariantsDFModerate_High_InOthers = SigBVariantsDFModerate_High %>% filter(BasePosition %in% Confounders_SigB$BasePosition)
SigBVariantsDFModerate_High_NotInOthers = SigBVariantsDFModerate_High %>% filter(!(BasePosition %in% Confounders_SigB$BasePosition))

Confounders_SigB = (Confounders_SigB %>% select(BasePosition, REF, ALT, Position, TotalPresent)) 
colnames(Confounders_SigB)= c("BasePosition", "RefUSA300", "AltOthers", "Position", "TotalPresent")
SigBVariantsDFModerate_High_InOthers = SigBVariantsDFModerate_High_InOthers %>% left_join(Confounders_SigB,by="BasePosition")

SigBVariantsDFModerate_High_InOthersSubset = SigBVariantsDFModerate_High_InOthers %>% select(Patient,TotalAltDepth, BasePosition, VariantType,  PredictedEffect, TimepointPresent, TimepointAbsent, RefUSA300, Ref,  AltOthers, Alt, Comparison, TotalPresent)


# Wherever there's a mismatch between USA300 genome reference alleles and the same base position of the , 
# it means that one alignment represents the reverse complement and the other represents the forward version of the gene
# (this is because I aligned fwd and reverse metagenomic reads separately)
# Correct for this by making RefCorrectedUSA300, AltCorrectedUSA300 columns which represent the complements of the snp-sites results,
# when there's a mismatch 
SigBVariantsDFModerate_High_InOthersSubset$RefCorrectedUSA300 <- sapply(1:nrow(SigBVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(SigBVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=SigBVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(SigBVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]), SigBVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"] ))
SigBVariantsDFModerate_High_InOthersSubset$AltCorrectedUSA300 <- sapply(1:nrow(SigBVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(SigBVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=SigBVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(SigBVariantsDFModerate_High_InOthersSubset[x,"AltOthers"]), SigBVariantsDFModerate_High_InOthersSubset[x,"AltOthers"] ))

# Now, make a column that represents Alt(the metagenomic SnpEff alt alleles listed for a position) but with any alts that were found in the 
# snp-sites results (which are from confounding species' sigB gene) removed; then, can filter out rows where only "<*>" remains. this will, effectively, 
# filter out any confounded variants (variants that might be due to these other species' presence in the metagenome)
SigBVariantsDFModerate_High_InOthersSubset$AltFiltered <- sapply(1:nrow(SigBVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) FilterConfounders(SigBVariantsDFModerate_High_InOthersSubset[x,"Alt"], SigBVariantsDFModerate_High_InOthersSubset[x,"AltCorrectedUSA300"]))
# Filter out any where only <*> is left after doing that. 
SigBVariantsDFModerate_High_InOthersSubsetFiltered = SigBVariantsDFModerate_High_InOthersSubset %>% filter(AltFiltered!="<*>")

SigBVariantsDFModerate_High_InOthersSubsetFiltered  = SigBVariantsDFModerate_High_InOthersSubsetFiltered %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

SigBVariantsDFModerate_High_NotInOthers  = SigBVariantsDFModerate_High_NotInOthers %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

SigBVariantsDFModerate_High_Filtered = rbind(SigBVariantsDFModerate_High_InOthersSubsetFiltered,SigBVariantsDFModerate_High_NotInOthers )

# Filter to range of SigB gene (no upstream variants please) and filter out indels 
SigBVariantsDFModerate_High_Filtered = SigBVariantsDFModerate_High_Filtered %>% filter(BasePosition %in% 2184996:2185766) %>% filter(!is.na(TotalAltDepth))


TestAltDepth <- function(fulldf, row){
  altdepth = row[2]
  altalleles = row[9]
  patient = row[1]
  position = row[3]
  present = row[6]
  comparison= row[10]
  altallelesList= str_split(altalleles, pattern=",")
  altallelesList = altallelesList[altallelesList!="<*>"]
  passes=F
  
  # If alt depth > # different alt alleles, set to true and return it
  if(altdepth > length(altallelesList)){
    passes=T
  }else{
    # otherwise, get all variants in that patient at that base position that aren't from the same 'timepoint present'
    smallDF = fulldf %>% filter( (BasePosition==position) & (Patient==patient) & (TimepointPresent!=present))
    if(nrow(smallDF) == 0){
      passes=F 
    }else{
      alleles_other = unlist(sapply(smallDF$Alt, function(x) str_split(x, pattern=",")))
      if(length(intersect(alleles_other,altallelesList )) > 0 ){
        passes=T
      }
    }
    
  }
  

  return(passes)
}




SigBVariantsDFModerate_High_Filtered$Passes= sapply(1:nrow(SigBVariantsDFModerate_High_Filtered), function(x) TestAltDepth(SigBVariantsDFModerate_High_Filtered, SigBVariantsDFModerate_High_Filtered[x,]))

FinalFilteredSigB = (SigBVariantsDFModerate_High_Filtered %>% filter(Passes==TRUE))
write.csv(FinalFilteredSigB, file="/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/FinalFilteredSigBVariants.csv")



#SigBVariantsDFModerate_High_Filtered = SigBVariantsDFModerate_High_Filtered %>% filter(BasePosition %in% 2184996:2185766) %>% filter(TotalAltDepth >=2) %>% arrange(Patient, BasePosition)

# Filter by allelic depth :

# 
# 
# 
# 
# 
# SigBVariantsDFModerate_High_Filtered_Gained = SigBVariantsDFModerate_High_Filtered %>% filter(TimepointPresent > TimepointAbsent)
# 
# SigBFilteredObs = data.frame(table(SigBVariantsDFModerate_High_Filtered_Gained[c("Patient","Comparison", "VariantType")]) )
# 
# SigBFilteredObs = SigBFilteredObs %>% filter(Freq>0) %>% arrange(Patient, Comparison)
# 
# SigBFilteredObs$EarlyTimepoint = sapply(SigBFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][1])
# SigBFilteredObs$LateTimepoint = sapply(SigBFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][2])
# 
# SigBFilteredObs %>% select(Patient, Comparison, EarlyTimepoint, LateTimepoint) %>% tidyr::expand(Patient, Comparison, EarlyTimepoint, LateTimepoint)

# SigBVariantsDFModerate_High_NotInOthers = SigBVariantsDFModerate_High_NotInOthers %>% filter(BasePosition %in% 2184996:2185766)
# SigBVariantsDFModerate_High_NotInOthers$EarlyTimepoint = sapply(SigBVariantsDFModerate_High_NotInOthers$Comparison, function(x) str_split(x, ":")[[1]][1])
# SigBVariantsDFModerate_High_NotInOthers$LateTimepoint = sapply(SigBVariantsDFModerate_High_NotInOthers$Comparison, function(x) str_split(x, ":")[[1]][2])
# SigBNotInOthers = SigBVariantsDFModerate_High_NotInOthers %>% select(Patient, Comparison, VariantType, EarlyTimepoint, )
# 

SigBFilteredObs

# RsbU
#########


Confounders_rsbU =read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/ConfoundersBlastN_Filtered/rsbUsnps", sep="\t", comment.char = "#", header=F)

Confounders_rsbU_old = Confounders_rsbU

Confounders_rsbU_old$BasePosition = 2187668-(Confounders_rsbU_old$V2 - 1)


Presence = Confounders_rsbU[, 10:ncol(Confounders_rsbU)]
Presence[Presence>0] <- 1
TotalPresent = rowSums(Presence)
Confounders_rsbU = Confounders_rsbU %>% select(V1,V2,V3,V4,V5,V6,V7,V8,V9)
Confounders_rsbU$TotalPresent = TotalPresent

colnames(Confounders_rsbU) = c("Chromosome", "Position", "ID","REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TotalPresent")
Confounders_rsbU$BasePosition=2187668-(Confounders_rsbU$Position - 1)

Comparisons_rsbU = Filter(function(x) grepl("rsbU",x), AdequateDepth)
rsbUcompareDF = data.frame(FileName=Comparisons_rsbU)

rsbUcompareDF$Patient = sapply(rsbUcompareDF$FileName, function(x) str_split(x, "_")[[1]][3])
unique(rsbUcompareDF$Patient)

rsbUVariantsDF = data.frame()
for(patient in unique(rsbUcompareDF$Patient)){
  subDF = rsbUcompareDF %>% filter(Patient==patient)
  
  for(i in 1:nrow(subDF)){
    fname=(subDF$FileName)[i]
    comparison= paste((str_split(fname, "_")[[1]][4]),(str_split(fname, "_")[[1]][5]),  sep=":")   
    
    ReadDF = read.table( paste0(comparisonpath, (subDF$FileName)[i]))
    ReadDF$Patient = patient
    ReadDF$Comparison = comparison
    rsbUVariantsDF = rbind(rsbUVariantsDF,ReadDF )
    
  }
  
}

rsbUVariantsDF$TotalAltDepth = sapply(rsbUVariantsDF$I16, function(x) sum(sapply(str_split(x, ",")[[1]][c(3,4)], as.numeric) ))



rsbUVariantsDFModerate_High = rsbUVariantsDF %>% filter(PredictedEffect!="LOW" )


rsbUVariantsDFModerate_High_InOthers = rsbUVariantsDFModerate_High %>% filter(BasePosition %in% Confounders_rsbU$BasePosition)
rsbUVariantsDFModerate_High_NotInOthers = rsbUVariantsDFModerate_High %>% filter(!(BasePosition %in% Confounders_rsbU$BasePosition))

Confounders_rsbU = (Confounders_rsbU %>% select(BasePosition, REF, ALT, Position, TotalPresent)) 
colnames(Confounders_rsbU)= c("BasePosition", "RefUSA300", "AltOthers", "Position", "TotalPresent")
rsbUVariantsDFModerate_High_InOthers = rsbUVariantsDFModerate_High_InOthers %>% left_join(Confounders_rsbU,by="BasePosition")

rsbUVariantsDFModerate_High_InOthersSubset = rsbUVariantsDFModerate_High_InOthers %>% select(Patient,TotalAltDepth, BasePosition, VariantType,  PredictedEffect, TimepointPresent, TimepointAbsent, RefUSA300, Ref,  AltOthers, Alt, Comparison, TotalPresent)


# Wherever there's a mismatch between USA300 genome reference alleles and the same base position of the , 
# it means that one alignment represents the reverse complement and the other represents the forward version of the gene
# (this is because I aligned fwd and reverse metagenomic reads separately)
# Correct for this by making RefCorrectedUSA300, AltCorrectedUSA300 columns which represent the complements of the snp-sites results,
# when there's a mismatch 
rsbUVariantsDFModerate_High_InOthersSubset$RefCorrectedUSA300 <- sapply(1:nrow(rsbUVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(rsbUVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=rsbUVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(rsbUVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]), rsbUVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"] ))
rsbUVariantsDFModerate_High_InOthersSubset$AltCorrectedUSA300 <- sapply(1:nrow(rsbUVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(rsbUVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=rsbUVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(rsbUVariantsDFModerate_High_InOthersSubset[x,"AltOthers"]), rsbUVariantsDFModerate_High_InOthersSubset[x,"AltOthers"] ))

# Now, make a column that represents Alt(the metagenomic SnpEff alt alleles listed for a position) but with any alts that were found in the 
# snp-sites results (which are from confounding species' rsbU gene) removed; then, can filter out rows where only "<*>" remains. this will, effectively, 
# filter out any confounded variants (variants that might be due to these other species' presence in the metagenome)
rsbUVariantsDFModerate_High_InOthersSubset$AltFiltered <- sapply(1:nrow(rsbUVariantsDFModerate_High_InOthersSubset),
                                                                 function(x) FilterConfounders(rsbUVariantsDFModerate_High_InOthersSubset[x,"Alt"], rsbUVariantsDFModerate_High_InOthersSubset[x,"AltCorrectedUSA300"]))
# Filter out any where only <*> is left after doing that. 
rsbUVariantsDFModerate_High_InOthersSubsetFiltered = rsbUVariantsDFModerate_High_InOthersSubset %>% filter(AltFiltered!="<*>")

rsbUVariantsDFModerate_High_InOthersSubsetFiltered  = rsbUVariantsDFModerate_High_InOthersSubsetFiltered %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

rsbUVariantsDFModerate_High_NotInOthers  = rsbUVariantsDFModerate_High_NotInOthers %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

rsbUVariantsDFModerate_High_Filtered = rbind(rsbUVariantsDFModerate_High_InOthersSubsetFiltered,rsbUVariantsDFModerate_High_NotInOthers )


rsbUVariantsDFModerate_High_Filtered = rsbUVariantsDFModerate_High_Filtered %>% filter(!is.na(TotalAltDepth))
rsbUVariantsDFModerate_High_Filtered$Passes= sapply(1:nrow(rsbUVariantsDFModerate_High_Filtered), function(x) TestAltDepth(rsbUVariantsDFModerate_High_Filtered, rsbUVariantsDFModerate_High_Filtered[x,]))

finalRsbUFiltered = rsbUVariantsDFModerate_High_Filtered %>% filter(Passes==T)

write.csv(finalRsbUFiltered, file="/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/FinalFilteredRsbUVariants.csv")





rsbUVariantsDFModerate_High_Filtered_Gained = rsbUVariantsDFModerate_High_Filtered %>% filter(TimepointPresent > TimepointAbsent)

rsbUVariantsDFModerate_High_Filtered %>% filter(TotalAltDepth >= 2)

rsbUFilteredObs = data.frame(table(rsbUVariantsDFModerate_High_Filtered_Gained[c("Patient","Comparison", "VariantType")]) )

rsbUFilteredObs = rsbUFilteredObs %>% filter(Freq>0) %>% arrange(Patient, Comparison)

rsbUFilteredObs$EarlyTimepoint = sapply(rsbUFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][1])
rsbUFilteredObs$LateTimepoint = sapply(rsbUFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][2])

rsbUFilteredObs %>% select(Patient, Comparison, EarlyTimepoint, LateTimepoint) %>% tidyr::expand(Patient, Comparison, EarlyTimepoint, LateTimepoint)


301*3


# rsbV
#########

Confounders_rsbV =read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/ConfoundersBlastN_Filtered/rsbVsnps", sep="\t", comment.char = "#", header=F)

Confounders_rsbV_old = Confounders_rsbV

Confounders_rsbV_old$BasePosition = 2186548-(Confounders_rsbV_old$V2 - 1)


Presence = Confounders_rsbV[, 10:ncol(Confounders_rsbV)]
Presence[Presence>0] <- 1
TotalPresent = rowSums(Presence)
Confounders_rsbV = Confounders_rsbV %>% select(V1,V2,V3,V4,V5,V6,V7,V8,V9)
Confounders_rsbV$TotalPresent = TotalPresent

colnames(Confounders_rsbV) = c("Chromosome", "Position", "ID","REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TotalPresent")
Confounders_rsbV$BasePosition=2186548-(Confounders_rsbV$Position - 1)


Comparisons_rsbV = Filter(function(x) grepl("rsbV",x), AdequateDepth)
rsbVcompareDF = data.frame(FileName=Comparisons_rsbV)
rsbVcompareDF$Patient = sapply(rsbVcompareDF$FileName, function(x) str_split(x, "_")[[1]][3])
unique(rsbVcompareDF$Patient)

rsbVVariantsDF = data.frame()
for(patient in unique(rsbVcompareDF$Patient)){
  subDF = rsbVcompareDF %>% filter(Patient==patient)
  for(i in 1:nrow(subDF)){
    fname=(subDF$FileName)[i]
    comparison= paste((str_split(fname, "_")[[1]][4]),(str_split(fname, "_")[[1]][5]),  sep=":")   
    
    ReadDF = read.table( paste0(comparisonpath, (subDF$FileName)[i]))
    ReadDF$Patient = patient
    ReadDF$Comparison=comparison
    rsbVVariantsDF = rbind(rsbVVariantsDF,ReadDF )
    
  }
  
}

rsbVVariantsDF$TotalAltDepth = sapply(rsbVVariantsDF$I16, function(x) sum(sapply(str_split(x, ",")[[1]][c(3,4)], as.numeric) ))

rsbVVariantsDFHigh = rsbVVariantsDF %>% filter(PredictedEffect=="HIGH" & TotalAltDepth>=2 )
rsbVVariantsDFModerate_High = rsbVVariantsDF %>% filter(PredictedEffect!="LOW" & TotalAltDepth>=2 )
rsbVVariantsDFModerate_High_InOthers = rsbVVariantsDFModerate_High %>% filter(BasePosition %in% Confounders_rsbV$BasePosition)

rsbVVariantsDFModerate_High_InOthers = rsbVVariantsDFModerate_High %>% filter(BasePosition %in% Confounders_rsbV$BasePosition)
rsbVVariantsDFModerate_High_NotInOthers = rsbVVariantsDFModerate_High %>% filter(!(BasePosition %in% Confounders_rsbV$BasePosition))

Confounders_rsbV = (Confounders_rsbV %>% select(BasePosition, REF, ALT, Position, TotalPresent)) 
colnames(Confounders_rsbV)= c("BasePosition", "RefUSA300", "AltOthers", "Position", "TotalPresent")
rsbVVariantsDFModerate_High_InOthers = rsbVVariantsDFModerate_High_InOthers %>% left_join(Confounders_rsbV,by="BasePosition")

rsbVVariantsDFModerate_High_InOthersSubset = rsbVVariantsDFModerate_High_InOthers %>% select(Patient,TotalAltDepth, BasePosition, VariantType,  PredictedEffect, TimepointPresent, TimepointAbsent, RefUSA300, Ref,  AltOthers, Alt, Comparison, TotalPresent)

# Wherever there's a mismatch between USA300 genome reference alleles and the same base position of the , 
# it means that one alignment represents the reverse complement and the other represents the forward version of the gene
# (this is because I aligned fwd and reverse metagenomic reads separately)
# Correct for this by making RefCorrectedUSA300, AltCorrectedUSA300 columns which represent the complements of the snp-sites results,
# when there's a mismatch 
rsbVVariantsDFModerate_High_InOthersSubset$RefCorrectedUSA300 <- sapply(1:nrow(rsbVVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(rsbVVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=rsbVVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(rsbVVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]), rsbVVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"] ))
rsbVVariantsDFModerate_High_InOthersSubset$AltCorrectedUSA300 <- sapply(1:nrow(rsbVVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(rsbVVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=rsbVVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(rsbVVariantsDFModerate_High_InOthersSubset[x,"AltOthers"]), rsbVVariantsDFModerate_High_InOthersSubset[x,"AltOthers"] ))

# Now, make a column that represents Alt(the metagenomic SnpEff alt alleles listed for a position) but with any alts that were found in the 
# snp-sites results (which are from confounding species' rsbV gene) removed; then, can filter out rows where only "<*>" remains. this will, effectively, 
# filter out any confounded variants (variants that might be due to these other species' presence in the metagenome)
rsbVVariantsDFModerate_High_InOthersSubset$AltFiltered <- sapply(1:nrow(rsbVVariantsDFModerate_High_InOthersSubset),
                                                                 function(x) FilterConfounders(rsbVVariantsDFModerate_High_InOthersSubset[x,"Alt"], rsbVVariantsDFModerate_High_InOthersSubset[x,"AltCorrectedUSA300"]))
# Filter out any where only <*> is left after doing that. 
rsbVVariantsDFModerate_High_InOthersSubsetFiltered = rsbVVariantsDFModerate_High_InOthersSubset %>% filter(AltFiltered!="<*>")

rsbVVariantsDFModerate_High_InOthersSubsetFiltered  = rsbVVariantsDFModerate_High_InOthersSubsetFiltered %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

rsbVVariantsDFModerate_High_NotInOthers  = rsbVVariantsDFModerate_High_NotInOthers %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

rsbVVariantsDFModerate_High_Filtered = rbind(rsbVVariantsDFModerate_High_InOthersSubsetFiltered,rsbVVariantsDFModerate_High_NotInOthers )

rsbVVariantsDFModerate_High_Filtered_Gained = rsbVVariantsDFModerate_High_Filtered %>% filter(TimepointPresent > TimepointAbsent)

rsbVFilteredObs = data.frame(table(rsbVVariantsDFModerate_High_Filtered_Gained[c("Patient","Comparison", "VariantType")]) )

rsbVFilteredObs = rsbVFilteredObs %>% filter(Freq>0) %>% arrange(Patient, Comparison)

rsbVFilteredObs$EarlyTimepoint = sapply(rsbVFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][1])
rsbVFilteredObs$LateTimepoint = sapply(rsbVFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][2])

rsbVFilteredObs %>% select(Patient, Comparison, EarlyTimepoint, LateTimepoint) %>% tidyr::expand(Patient, Comparison, EarlyTimepoint, LateTimepoint)

rsbVVariantsDFModerate_High_Filtered = rsbVVariantsDFModerate_High_Filtered %>% filter(!is.na(TotalAltDepth))
rsbVVariantsDFModerate_High_Filtered$Passes= sapply(1:nrow(rsbVVariantsDFModerate_High_Filtered), function(x) TestAltDepth(rsbVVariantsDFModerate_High_Filtered, rsbVVariantsDFModerate_High_Filtered[x,]))

finalRsbVFiltered = rsbVVariantsDFModerate_High_Filtered %>% filter(Passes==T)

write.csv(finalRsbVFiltered, file="/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/FinalFilteredRsbVVariants.csv")





# rsbW
#########

Confounders_rsbW =read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/ConfoundersBlastN_Filtered/rsbWsnps", sep="\t", comment.char = "#", header=F)

Confounders_rsbW_old = Confounders_rsbW

Confounders_rsbW_old$BasePosition = 2186220-(Confounders_rsbW_old$V2 - 1)


Presence = Confounders_rsbW[, 10:ncol(Confounders_rsbW)]
Presence[Presence>0] <- 1
TotalPresent = rowSums(Presence)
Confounders_rsbW = Confounders_rsbW %>% select(V1,V2,V3,V4,V5,V6,V7,V8,V9)
Confounders_rsbW$TotalPresent = TotalPresent

colnames(Confounders_rsbW) = c("Chromosome", "Position", "ID","REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TotalPresent")
Confounders_rsbW$BasePosition=2186220-(Confounders_rsbW$Position - 1)


Comparisons_rsbW = Filter(function(x) grepl("rsbW",x), AdequateDepth)
rsbWcompareDF = data.frame(FileName=Comparisons_rsbW)
rsbWcompareDF$Patient = sapply(rsbWcompareDF$FileName, function(x) str_split(x, "_")[[1]][3])
unique(rsbWcompareDF$Patient)

rsbWVariantsDF = data.frame()
for(patient in unique(rsbWcompareDF$Patient)){
  subDF = rsbWcompareDF %>% filter(Patient==patient)
  for(i in 1:nrow(subDF)){
    fname=(subDF$FileName)[i]
    comparison= paste((str_split(fname, "_")[[1]][4]),(str_split(fname, "_")[[1]][5]),  sep=":")   
    
    ReadDF = read.table( paste0(comparisonpath, (subDF$FileName)[i]))
    ReadDF$Patient = patient
    ReadDF$Comparison=comparison
    rsbWVariantsDF = rbind(rsbWVariantsDF,ReadDF )
    
  }
  
}

rsbWVariantsDF$TotalAltDepth = sapply(rsbWVariantsDF$I16, function(x) sum(sapply(str_split(x, ",")[[1]][c(3,4)], as.numeric) ))

rsbWVariantsDFHigh = rsbWVariantsDF %>% filter(PredictedEffect=="HIGH" & TotalAltDepth>=2 )
rsbWVariantsDFModerate_High = rsbWVariantsDF %>% filter(PredictedEffect!="LOW" & TotalAltDepth>=2 )
rsbWVariantsDFModerate_High_InOthers = rsbWVariantsDFModerate_High %>% filter(BasePosition %in% Confounders_rsbW$BasePosition)

rsbWVariantsDFModerate_High_InOthers = rsbWVariantsDFModerate_High %>% filter(BasePosition %in% Confounders_rsbW$BasePosition)
rsbWVariantsDFModerate_High_NotInOthers = rsbWVariantsDFModerate_High %>% filter(!(BasePosition %in% Confounders_rsbW$BasePosition))

Confounders_rsbW = (Confounders_rsbW %>% select(BasePosition, REF, ALT, Position, TotalPresent)) 
colnames(Confounders_rsbW)= c("BasePosition", "RefUSA300", "AltOthers", "Position", "TotalPresent")
rsbWVariantsDFModerate_High_InOthers = rsbWVariantsDFModerate_High_InOthers %>% left_join(Confounders_rsbW,by="BasePosition")

rsbWVariantsDFModerate_High_InOthersSubset = rsbWVariantsDFModerate_High_InOthers %>% select(Patient,TotalAltDepth, BasePosition, VariantType,  PredictedEffect, TimepointPresent, TimepointAbsent, RefUSA300, Ref,  AltOthers, Alt, Comparison, TotalPresent)

# Wherever there's a mismatch between USA300 genome reference alleles and the same base position of the , 
# it means that one alignment represents the reverse complement and the other represents the forward version of the gene
# (this is because I aligned fwd and reverse metagenomic reads separately)
# Correct for this by making RefCorrectedUSA300, AltCorrectedUSA300 columns which represent the complements of the snp-sites results,
# when there's a mismatch 
rsbWVariantsDFModerate_High_InOthersSubset$RefCorrectedUSA300 <- sapply(1:nrow(rsbWVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(rsbWVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=rsbWVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(rsbWVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]), rsbWVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"] ))
rsbWVariantsDFModerate_High_InOthersSubset$AltCorrectedUSA300 <- sapply(1:nrow(rsbWVariantsDFModerate_High_InOthersSubset),
                                                                        function(x) if_else(rsbWVariantsDFModerate_High_InOthersSubset[x,"RefUSA300"]!=rsbWVariantsDFModerate_High_InOthersSubset[x,"Ref"],
                                                                                            GetComplement(rsbWVariantsDFModerate_High_InOthersSubset[x,"AltOthers"]), rsbWVariantsDFModerate_High_InOthersSubset[x,"AltOthers"] ))

# Now, make a column that represents Alt(the metagenomic SnpEff alt alleles listed for a position) but with any alts that were found in the 
# snp-sites results (which are from confounding species' rsbW gene) removed; then, can filter out rows where only "<*>" remains. this will, effectively, 
# filter out any confounded variants (variants that might be due to these other species' presence in the metagenome)
rsbWVariantsDFModerate_High_InOthersSubset$AltFiltered <- sapply(1:nrow(rsbWVariantsDFModerate_High_InOthersSubset),
                                                                 function(x) FilterConfounders(rsbWVariantsDFModerate_High_InOthersSubset[x,"Alt"], rsbWVariantsDFModerate_High_InOthersSubset[x,"AltCorrectedUSA300"]))
# Filter out any where only <*> is left after doing that. 
rsbWVariantsDFModerate_High_InOthersSubsetFiltered = rsbWVariantsDFModerate_High_InOthersSubset %>% filter(AltFiltered!="<*>")

rsbWVariantsDFModerate_High_InOthersSubsetFiltered  = rsbWVariantsDFModerate_High_InOthersSubsetFiltered %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

rsbWVariantsDFModerate_High_NotInOthers  = rsbWVariantsDFModerate_High_NotInOthers %>% select(Patient, TotalAltDepth, BasePosition, VariantType, PredictedEffect, TimepointPresent, TimepointAbsent, Ref, Alt, Comparison)

rsbWVariantsDFModerate_High_Filtered = rbind(rsbWVariantsDFModerate_High_InOthersSubsetFiltered,rsbWVariantsDFModerate_High_NotInOthers )

rsbWVariantsDFModerate_High_Filtered_Gained = rsbWVariantsDFModerate_High_Filtered %>% filter(TimepointPresent > TimepointAbsent)

rsbWFilteredObs = data.frame(table(rsbWVariantsDFModerate_High_Filtered_Gained[c("Patient","Comparison", "VariantType")]) )

rsbWFilteredObs = rsbWFilteredObs %>% filter(Freq>0) %>% arrange(Patient, Comparison)

rsbWFilteredObs$EarlyTimepoint = sapply(rsbWFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][1])
rsbWFilteredObs$LateTimepoint = sapply(rsbWFilteredObs$Comparison, function(x) str_split(x, ":")[[1]][2])

rsbWFilteredObs %>% select(Patient, Comparison, EarlyTimepoint, LateTimepoint) %>% tidyr::expand(Patient, Comparison, EarlyTimepoint, LateTimepoint)



rsbWVariantsDFModerate_High_Filtered = rsbWVariantsDFModerate_High_Filtered %>% filter(!is.na(TotalAltDepth))
rsbWVariantsDFModerate_High_Filtered$Passes= sapply(1:nrow(rsbWVariantsDFModerate_High_Filtered), function(x) TestAltDepth(rsbWVariantsDFModerate_High_Filtered, rsbWVariantsDFModerate_High_Filtered[x,]))

finalRsbWFiltered = rsbWVariantsDFModerate_High_Filtered %>% filter(Passes==T)

write.csv(finalRsbWFiltered, file="/Users/amycampbell/Documents/DataInputGithub/data/VariantsSigB/FinalFilteredRsbWVariants.csv")

