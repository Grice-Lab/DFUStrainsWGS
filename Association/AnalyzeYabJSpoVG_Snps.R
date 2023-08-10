# Amy Campbell
# 2023
# aligning yabJ-spoVG operon genes gene-by-gene
# analyzing snp-sites output of that alignment

library(Biostrings) # Biostrings_2.62.0
library(dplyr) # dplyr_1.1.1
library(stringr) # stringr_1.5.0
library(ggplot2) # 3.4.2 


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


