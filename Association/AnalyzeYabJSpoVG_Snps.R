# Amy Campbell
# 2023
# aligning yabJ-spoVG operon genes gene-by-gene
# analyzing snp-sites output of that alignment

library(Biostrings) # Biostrings_2.62.0
library(dplyr) # dplyr_1.1.1
library(stringr) # stringr_1.5.0
library(ggplot2) # 3.4.2 

Phenotypes = read.csv("~/Desktop/GriceLabGit/Staphyloxanthin/Data/InVitroData/staphyloxanthin_paper_updated.csv")
Phenotypes$DORN = sapply(Phenotypes$IsolateID, function(x) str_replace(x,"SA", "DORN"))
CCmap = read.csv("Documents/DataInputGithub/data/Phylogeny2022Data/CCMapPlotting.csv")
Patient_Genome_Info = read.csv("~/Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")
Patient_Genome_Info$DORN = paste0("DORN", Patient_Genome_Info$Doern.lab.bank.)
patient_genome = Patient_Genome_Info %>% select(patient_id, DORN)

yabJsnps = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/yabJ_snps.txt"
yabJsequence = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/yabJORF_fpr3757.fasta"
yabJalignment = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/yabJ_mafft.aln"

spoVGsnps = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/spoVG_snps.txt"
spoVGsequence = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/spoVG_fpr3757.fasta"
spoVGalignment = "/Users/amycampbell/Documents/DataInputGithub/data/yabJspoVG/spoVG_mafft.aln"

######################
# yabJ alignment SNPs
######################
# There's one deletion snp (just in SA1086) but 18 different positions with base substitutions
# There's also one position where there's two different alt alleles 

YabJresults = read.table(yabJsnps,skip=3,sep='\t', comment.char="$",header=T)
refSequenceYabJ = Biostrings::readDNAStringSet(yabJsequence, format="fasta")
alignmentSequencesYabJ=(Biostrings::readBStringSet(yabJalignment, format="fasta"))$yabJ_fpr3757

alignmentbases = str_split(as.character(unlist(alignmentSequencesYabJ)), "")[[1]]
alignmentbases = sapply(alignmentbases, toupper)

indexer =1 
alignmentpositions = 1:length(alignmentbases)
refpositions = c()

#Iterate through the bases of the alignment
for(b in alignmentbases ){
  # If the alignment version of the ref genome="*", then there's a gap in it there. 
  # Append "NA" to the ref positions since the base in the reference genome at that point is NA
  # this deals with things that are missing from the reference genome (and are therefore gaps)
  if(b == "*"){
    refpositions = append(refpositions, "NA")
  }else{
    refpositions = append(refpositions, indexer)
    indexer=indexer+1
  }
}


RefBasesYabJ = str_split(as.character(unlist(refSequenceYabJ)), "")[[1]]

AANums =sapply(1:length(RefBasesYabJ), function(x) ceiling(x/3))
refDF = data.frame(RefPos = 1:length(RefBasesYabJ), refbase = RefBasesYabJ, AAnum =AANums, codonPos = rep(c(1,2,3), length(RefBasesYabJ)/3))
refDF$RefPos = sapply(refDF$RefPos, as.character)
MapRef_Alignment = data.frame(POS = alignmentpositions,RefPos = refpositions)

MapRef_Alignment$RefPos = sapply(MapRef_Alignment$RefPos, as.character)

refDF = MapRef_Alignment %>% left_join(refDF,by="RefPos")
row.names(refDF) = refDF$POS


newResults = data.frame()
for(j in 1:nrow(YabJresults)){
  rowitem = YabJresults[j,]
  
  if(grepl(pattern=",", rowitem["ALT"])){
    altcharslist = str_split(rowitem["ALT"], ",")[[1]]

    for(ind in 1:length(altcharslist)){
      newrow=rowitem
      newrow_numeric = newrow[10:length(newrow)]
      newrow_numeric[newrow_numeric!=ind] = 0
      newrow_numeric[newrow_numeric==ind] = 1
      newrow[10:length(newrow)] <- newrow_numeric
      newrow["ALT"] <- altcharslist[ind]
      newResults=rbind(newResults, newrow)
    }

  }else{
    newResults = rbind(newResults,rowitem)
    
  }
}


typelist=c()
describelist=c()

for(i in 1:nrow(newResults)){
  SNPinfo =newResults[i,]
  # We know there's only one deletion mutation and it's in the alt, not the ref
  # We also know it's a frameshift
  # If we were dealing with a different alignment where there were gaps in the reference genome,
  # Would have to add more here
  if(SNPinfo["ALT"]=="*"){
    type_mute ="FS"
    description="FrameshiftMutation"
  }else{
    refcodonrows = refDF %>% filter(AAnum==refDF[as.character(SNPinfo["POS"]),"AAnum"])
    row.names(refcodonrows) = refcodonrows$POS
    refcodon_seq = paste(refcodonrows$refbase, collapse="")

    refcodonrows$AltBase = NA
    refcodonrows[as.character(SNPinfo["POS"]), "AltBase"] <- SNPinfo["ALT"]
    refcodonrows$AltBase = if_else(is.na(refcodonrows$AltBase), refcodonrows$refbase, refcodonrows$AltBase)
    altcodon_seq = paste(refcodonrows$AltBase, collapse="")
    refAA =  Biostrings::GENETIC_CODE[[refcodon_seq]]
    altAA =  Biostrings::GENETIC_CODE[[altcodon_seq]]
    if(refAA==altAA){
      type_mute="Synonymous"
      description="SynonymousMutation"
    }else{
      if(altAA=="*"){
        type_mute="Stop_Gained"
        description=paste(refAA, "to", "STOP", sep="_")
        
      } else{
        type_mute="AASub"
        description=paste(refAA, "to", altAA, sep="_")
      }
    }
  }
  typelist = append(typelist, type_mute)
  describelist = append(describelist, description)
}

newResults$Type = typelist
newResults$Description = describelist


######################
# spoVG
######################

spoVGresults = read.table(spoVGsnps,skip=3,sep='\t', comment.char="$",header=T)
spoVGresults$ALT[spoVGresults$ALT==TRUE] <- "T"
refSequencespoVG= Biostrings::readDNAStringSet(spoVGsequence, format="fasta")
alignmentSequencesspoVG=(Biostrings::readBStringSet(spoVGalignment, format="fasta"))$spoVG_fpr3757

alignmentbases = str_split(as.character(unlist(alignmentSequencesspoVG)), "")[[1]]
alignmentbases = sapply(alignmentbases, toupper)

indexer =1 
alignmentpositions = 1:length(alignmentbases)
refpositions = c()

#Iterate through the bases of the alignment
for(b in alignmentbases ){
  # If the alignment version of the ref genome="*", then there's a gap in it there. 
  # Append "NA" to the ref positions since the base in the reference genome at that point is NA
  # this deals with things that are missing from the reference genome (and are therefore gaps)
  if(b == "*"){
    refpositions = append(refpositions, "NA")
  }else{
    refpositions = append(refpositions, indexer)
    indexer=indexer+1
  }
}


RefBasesSpoVG = str_split(as.character(unlist(refSequencespoVG)), "")[[1]]

AANums =sapply(1:length(RefBasesSpoVG), function(x) ceiling(x/3))
refDF = data.frame(RefPos = 1:length(RefBasesSpoVG), refbase = RefBasesSpoVG, AAnum =AANums, codonPos = rep(c(1,2,3), length(RefBasesSpoVG)/3))
refDF$RefPos = sapply(refDF$RefPos, as.character)
MapRef_Alignment = data.frame(POS = alignmentpositions,RefPos = refpositions)

MapRef_Alignment$RefPos = sapply(MapRef_Alignment$RefPos, as.character)

refDF = MapRef_Alignment %>% left_join(refDF,by="RefPos")
row.names(refDF) = refDF$POS


newResultsSpoVG = data.frame()
for(j in 1:nrow(spoVGresults)){
  rowitem = spoVGresults[j,]
  if(grepl(pattern=",", rowitem["ALT"])){
    altcharslist = str_split(rowitem["ALT"], ",")[[1]]
    for(ind in 1:length(altcharslist)){
      newrow=rowitem
      newrow_numeric = newrow[10:length(newrow)]
      newrow_numeric[newrow_numeric!=ind] = 0
      newrow_numeric[newrow_numeric==ind] = 1
      newrow[10:length(newrow)] <- newrow_numeric
      newrow["ALT"] <- altcharslist[ind]
      newResultsSpoVG=rbind(newResultsSpoVG, newrow)
    }
    
  }else{
    newResultsSpoVG = rbind(newResultsSpoVG,rowitem)
  }
}


typelist=c()
describelist=c()

for(i in 1:nrow(newResultsSpoVG)){
  SNPinfo =newResultsSpoVG[i,]
  # We know there's only one deletion mutation and it's in the alt, not the ref
  # We also know it's a frameshift
  # If we were dealing with a different alignment where there were gaps in the reference genome,
  # Would have to add more here
  if(SNPinfo["ALT"]=="*"){
    type_mute ="FS"
    description="FrameshiftMutation"
  }else{
    refcodonrows = refDF %>% filter(AAnum==refDF[as.character(SNPinfo["POS"]),"AAnum"])
    row.names(refcodonrows) = refcodonrows$POS
    refcodon_seq = paste(refcodonrows$refbase, collapse="")
    
    refcodonrows$AltBase = NA
    refcodonrows[as.character(SNPinfo["POS"]), "AltBase"] <- SNPinfo["ALT"]
    refcodonrows$AltBase = if_else(is.na(refcodonrows$AltBase), refcodonrows$refbase, refcodonrows$AltBase)
    altcodon_seq = paste(refcodonrows$AltBase, collapse="")
    refAA =  Biostrings::GENETIC_CODE[[refcodon_seq]]
    altAA =  Biostrings::GENETIC_CODE[[altcodon_seq]]
    if(refAA==altAA){
      type_mute="Synonymous"
      description="SynonymousMutation"
    }else{
      if(altAA=="*"){
        type_mute="Stop_Gained"
        description=paste(refAA, "to", "STOP", sep="_")
        
      } else{
        type_mute="AASub"
        description=paste(refAA, "to", altAA, sep="_")
      }
    }
  }
  typelist = append(typelist, type_mute)
  describelist = append(describelist, description)
}

newResultsSpoVG$Type = typelist
newResultsSpoVG$Description = describelist
newResultsSpoVG$Gene = "SpoVG"
newResults$Gene = "YabJ"


newResultsSpoVG$spoVG_fpr3757 = NULL
newResults$yabJ_fpr3757 = NULL

refDF

# Issue here:
# DORN925  is positive for an AA subtitution at POS115 corresponding to H-->D
# However, when you actually do a blastP, this position has an E (so GAA or GAG as opposed to GAT/GAC)


yabJspoVGORFs =rbind(newResults, newResultsSpoVG)
yabJspoVGORFs$gene_pos_desc = paste(yabJspoVGORFs$Gene, yabJspoVGORFs$POS, yabJspoVGORFs$Description, sep="_")

VariantPresence_Absence = yabJspoVGORFs %>% select(gene_pos_desc ,colnames(yabJspoVGORFs)[grepl("DORN", colnames(yabJspoVGORFs))])

NS_Only = VariantPresence_Absence %>% filter(!grepl("Synonymous", gene_pos_desc))

newcols = NS_Only$gene_pos_desc
NS_Only$gene_pos_desc=NULL
Summary_By_Genome= data.frame(t(NS_Only))
colnames(Summary_By_Genome) = newcols
Summary_By_Genome$DORN = rownames(Summary_By_Genome)

Summary_By_Genome = Summary_By_Genome %>% left_join(Phenotypes %>% select(DORN, staphyloxanthin, biofilm),by="DORN")
Summary_By_Genome = Summary_By_Genome %>% left_join(CCmap, by="DORN")
Summary_By_Genome = Summary_By_Genome %>% left_join(Patient_Genome_Info %>% select(patient_id, DORN), by="DORN")


Summary_By_Genome$SumVariants = rowSums(Summary_By_Genome[1:(ncol(Summary_By_Genome)-5)])
Summary_By_Genome$AnyVariants = if_else(Summary_By_Genome$SumVariants > 0, "Yes", "No")


# First, figure out how many different patients a variant was found in
# Then, figure out how many different CCs it was in 


Summary_By_Genome$CCpatient =paste(Summary_By_Genome$CCLabel, Summary_By_Genome$patient_id, sep="_")
results_by_variant = data.frame()


# look for variants which are variably present within CC, patient combinations 
VariableWithinCC_patient = data.frame()

# Look for variants which are variably present within CC (though maybe not within patinet)
VariableWithinCC = data.frame()

for(colid in newcols){
  Variant = Summary_By_Genome[, c(colid,"DORN", "staphyloxanthin", "biofilm","CCLabel", "patient_id", "CCpatient")]
  colnames(Variant) = c("VariantPresence", "DORN", "staphyloxanthin", "biofilm","CCLabel", "patient_id", "CCpatient")
  VariantPositive= Variant %>% filter(VariantPresence==1)
  VariantNegative= Variant %>% filter(VariantPresence==0)
  
  numgenomes = nrow(VariantPositive)
  NumCCs = length(unique(VariantPositive$CCLabel))
  NumPatients =  length(unique(VariantPositive$patient_id))
  CCs = paste(unique(VariantPositive$CCLabel), collapse=";")
  VariantPositive$CCpatient = paste(VariantPositive$CCLabel, VariantPositive$patient_id, sep="_")
  meanSTXpositive = mean(VariantPositive$staphyloxanthin, na.rm=T)
  meanSTXnegative = mean(VariantNegative$staphyloxanthin, na.rm=T)
  
  meanBiofilmpositive=mean(VariantPositive$biofilm, na.rm=T)
  meanBiofilmnegative=mean(VariantNegative$biofilm, na.rm=T)
  
  results_by_variant = rbind(results_by_variant, c(colid,NumPatients,NumCCs, CCs, meanSTXpositive, meanSTXnegative,meanBiofilmpositive, meanBiofilmnegative  ))
  for(ccpatientcombo in unique(VariantPositive$CCpatient)){
    NumPosiPatient = nrow(VariantPositive %>% filter(CCpatient == ccpatientcombo))
    NumNegativePatient = nrow(VariantNegative %>% filter(CCpatient == ccpatientcombo))
    if(NumPosiPatient>0 & NumNegativePatient > 0){
      VariableWithinCC_patient = rbind(VariableWithinCC_patient, c(colid,NumPosiPatient, NumNegativePatient, ccpatientcombo ))
    }
  }
  for(CC in unique(VariantPositive$CCLabel)){
    NumPosi = nrow(VariantPositive %>% filter(CCLabel == CC))
    numpatientsposi = length(unique( (VariantPositive %>% filter(CCLabel == CC))$patient_id))
    numpatientsnegative = length(unique( (VariantNegative %>% filter(CCLabel == CC))$patient_id))
    
    NumNegative = nrow(VariantNegative %>% filter(CCLabel == CC))
    if((NumPosi>0) & (NumNegative > 0)){
      VariableWithinCC = rbind(VariableWithinCC, c(colid,NumPosi, NumNegative, CC, numpatientsposi,numpatientsnegative  ))
    }
  }
}


colnames(results_by_variant) = c("Variant", "Number Patients", "Number CCs", "CCs", "meanSTXpositive", "meanSTXnegative", "meanBiofilmpositive", "meanBiofilmnegative")
colnames(VariableWithinCC) = c("Variant", "NumberGenomesWithVariant_CC", "NumberGenomesWithout_CC","CC", "NumberPatientsWith", "NumberPatientsWithout")




# The multiallelic variant at AA117
###################################
AA117alleles = Summary_By_Genome %>% select(YabJ_117_H_to_Q, staphyloxanthin, CCLabel)

ggplot(AA117alleles, aes(x=factor(YabJ_117_H_to_Q), y=log(staphyloxanthin))) + geom_boxplot(fill="#B8860B") + geom_jitter()


# YabJ 115 (H-->D)
# YabJ 117 (H-->Q)
# (the two are identically distributed)
#######################################
YabJ_CC1s_115_117_HtoDQ = Summary_By_Genome %>% select(YabJ_115_H_to_D, YabJ_117_H_to_Q, CCLabel, patient_id, staphyloxanthin, biofilm) %>% filter(CCLabel == "CC1")
all(YabJ_CC1s_115_117_HtoDQ$YabJ_115_H_to_D==YabJ_CC1s_115_117_HtoDQ$YabJ_117_H_to_Q)

ggplot(YabJ_CC1s_115_117_HtoDQ, aes(x=factor(YabJ_117_H_to_Q),y=log10(staphyloxanthin))) + geom_boxplot() + geom_jitter() 
ggplot(YabJ_CC1s_115_117_HtoDQ, aes(x=factor(YabJ_117_H_to_Q),y=log10(biofilm))) + geom_boxplot() + geom_jitter() 




YabJ_251 = Summary_By_Genome %>% select(YabJ_251_FrameshiftMutation, CCLabel, patient_id, staphyloxanthin, biofilm) %>% filter(CCLabel == "CC1" & patient_id=="310")
ggplot(YabJ_251, aes(x=factor(YabJ_251_FrameshiftMutation),y=log10(staphyloxanthin))) + geom_boxplot() + geom_jitter() 
ggplot(YabJ_251, aes(x=factor(YabJ_251_FrameshiftMutation),y=log10(biofilm))) + geom_boxplot() + geom_jitter() 

CC5_SpoVG_160_D_to_Y = Summary_By_Genome %>% select(SpoVG_160_D_to_Y, CCLabel, patient_id, staphyloxanthin, biofilm) %>% filter(CCLabel == "CC5")
ggplot(CC5_SpoVG_160_D_to_Y, aes(x=factor(SpoVG_160_D_to_Y),y=log10(staphyloxanthin))) + geom_boxplot() + geom_jitter() 
ggplot(CC5_SpoVG_160_D_to_Y, aes(x=factor(SpoVG_160_D_to_Y),y=log10(biofilm))) + geom_boxplot() + geom_jitter() 

CC5_SpoVG_290_S_to_L = Summary_By_Genome %>% select(SpoVG_290_S_to_L, CCLabel, patient_id, staphyloxanthin, biofilm) %>% filter(CCLabel == "CC5")
ggplot(CC5_SpoVG_290_S_to_L, aes(x=factor(SpoVG_290_S_to_L),y=log10(staphyloxanthin))) + geom_boxplot() + geom_jitter() 
ggplot(CC5_SpoVG_290_S_to_L, aes(x=factor(SpoVG_290_S_to_L),y=log10(biofilm))) + geom_boxplot() + geom_jitter() 






# Identifying consistent variants between 1086 and (1082 + 1081)
#################################################################

# NucDiff comparisons of 1086 (as the ref) to 1081 
##################################################
DORN1081 = read.csv2("Documents/DataInputGithub/data/yabJspoVG/nucdiff_patient310/DORN1086_DORN1081/results/DORN1086_DORN1081_ref_snps.gff" ,comment.char = "#", sep="\t",header=F)
colnames(DORN1081) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
DORN1081 = DORN1081 %>% select(Start, End, Description)

Structural1081 = read.csv2("Documents/DataInputGithub/data/yabJspoVG/nucdiff_patient310/DORN1086_DORN1081/results/DORN1086_DORN1081_ref_struct.gff" ,comment.char = "#", sep="\t",header=F)
colnames(Structural1081) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Structural1081 = Structural1081 %>% select(Start, End, Description)

Additional1081 = read.csv2("Documents/DataInputGithub/data/yabJspoVG/nucdiff_patient310/DORN1086_DORN1081/results/DORN1086_DORN1081_ref_additional.gff" ,comment.char = "#", sep="\t",header=F)
colnames(Additional1081) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Additional1081 = Additional1081 %>% select(Start, End, Description)

allvariants1081 = rbind(DORN1081,Structural1081, Additional1081)



# NucDiff comparisons of 1086 (ref) to 1082 
############################################
DORN1082 = read.csv2("Documents/DataInputGithub/data/yabJspoVG/nucdiff_patient310/DORN1086_DORN1082/results/DORN1086_DORN1082_ref_snps.gff" ,comment.char = "#", sep="\t",header=F)
colnames(DORN1082) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
DORN1082 = DORN1082 %>% select(Start, End, Description)

Structural1082 = read.csv2("Documents/DataInputGithub/data/yabJspoVG/nucdiff_patient310/DORN1086_DORN1082/results/DORN1086_DORN1082_ref_struct.gff" ,comment.char = "#", sep="\t",header=F)
colnames(Structural1082) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Structural1082 = Structural1082 %>% select(Start, End, Description)


Additional1082 = read.csv2("Documents/DataInputGithub/data/yabJspoVG/nucdiff_patient310/DORN1086_DORN1082/results/DORN1086_DORN1082_ref_additional.gff" ,comment.char = "#", sep="\t",header=F)
colnames(Additional1082) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Additional1082 = Additional1082 %>% select(Start, End, Description)
allvariants1082 = rbind(DORN1082,Structural1082, Additional1082)


allvariants1082$StartStop = paste(allvariants1082$Start, allvariants1082$End, sep="_")
allvariants1081$StartStop = paste(allvariants1081$Start, allvariants1081$End, sep="_")

View(allvariants1081 %>% filter(StartStop %in% (intersect(allvariants1082$StartStop, allvariants1081$StartStop))))

