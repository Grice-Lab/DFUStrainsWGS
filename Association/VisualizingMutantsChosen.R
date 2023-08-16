# Amy Campbell
# KDE subset GWAS results
# Based on code from spring 2022 -- updated to just include used code in 2023

library(tidyr)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(ggpubr)


################
# FUNCTIONS
################
FindAnnotation <- function(unitigposition,unitiglength, AnnotationFrame,ContigName,unitigstrand, returnType){
  # INPUTS:
  # - unitigposition - starting position of unitig 
  # - unitiglength - sequence length of the unitig in question
  # - AnnotationFrame - the GFF data frame of the annotated reference genome 
  # - returnType
  # 1 : Annotation value; tells the annotation value
  # 2 : Overlaps ; tells how many bp this unitig overlaps with an actual CDS annotation by
  # 3:  Tells the start locus of the closest annotation (useful in cases of no overlap)
  # 4: Tells the number of genes it overlapped with (like in the case of a unitig that is huge and represents a whole insertion or something)
  # OUTPUTS:
  # Either full annotation string from .gff located at the CDS start closest to the unitig start
  # or the # of overlapping basepairs between the annotation and unitig depending on returnType val
  
  # Unitig start
  unitigposition = (as.integer(unitigposition))
  
  # Unitig length
  unitiglength = as.integer(unitiglength)
  
  # Otherwise, look for a row where any part of the start:end interval in that row 
  # overlaps with any part of the unitig's start:end interval 
  if(unitigstrand =="-"){
    UnitigStart = unitigposition - unitiglength + 1 
    UnitigEnd = unitigposition
  }else if(unitigstrand=="+"){
    UnitigStart = unitigposition
    UnitigEnd = unitigposition + unitiglength - 1 
  }
  
  
  # Check for gene with start closest to the end of the unitig 
  closest_start_Unitig_start = which.min( abs(AnnotationFrame$start-UnitigEnd) )
  RowClosest = AnnotationFrame[closest_start_Unitig_start,]
  
  RowsPossible = AnnotationFrame %>% filter(seqname==ContigName)
  
  # Overlaps are possible wherever BOTH the Unitig starts before the end of the gene
  # AND the unitig ends after or at the start of the gene 
  RowsPossible = RowsPossible %>% filter((UnitigStart <= end) & (UnitigEnd >= start))
  
  # Return either the overlapping annotation or the annotation the end of the unitig is closest to the start of
  if(returnType == 1){
    if(nrow(RowsPossible)==0){
      
      return(RowClosest["attr"])
      
    } else if(nrow(RowsPossible)==1){
      
      return(RowsPossible["attr"])
    } else if(nrow(RowsPossible) >1){
      return(RowsPossible[1, "attr"])
    }
  }
  else if(returnType==2){
    if(nrow(RowsPossible)==0){
      return(0)
    }else if(nrow(RowsPossible)==1){
      
      # If the unitig starts before the start of the annotation (e.g. if UnitigStart < start ), the overlap will be genestart:unitig end  
      # If the unitig starts somewhere in the middle of the annotation (e.g. UnitigStart > start), the overlap will be min(unitigend - unitigstart, end-unitigstart)
      if(UnitigStart < RowsPossible[1,"start"] ){
        return( (UnitigEnd - RowsPossible[1,"start"]) )
      }else{
        
        return( min( (RowsPossible[1,"end"] - UnitigStart ), unitiglength) )
      }
      
    } else if(nrow(RowsPossible)>1){
      TotalRegionStart = RowsPossible[1, "start"]
      TotalRegionEnd = RowsPossible[nrow(RowsPossible),"end"]
      if(UnitigStart < TotalRegionStart){
        return(UnitigEnd - TotalRegionStart)
      }else{
        Result = (min(TotalRegionEnd - UnitigStart, unitiglength))
        return(Result[1])
      }
    }
    
  }
  else if(returnType==3){
    if(nrow(RowsPossible)==0){
      return(RowClosest[1,"start"])
      
    }else{
      return(RowsPossible[1, "start"])
    }
    
  }
  else if(returnType==4){
    return(nrow(RowsPossible))
  }
  
}





# Visualizing phenotypic relationships of mutants chosen for follow-up

ColorsCC = read.csv("~/Documents/DataInputGithub/data/Phylogeny2022Data/CC_Colorscheme.csv")

setwd("~/Documents/DataInputGithub/")

CCmap = read.csv("data/Phylogeny2022Data/CCMapPlotting.csv")
DFUIsolatesInfo = read.csv("data/DFU_Staph_aureus_isolates.csv")
DFUIsolatesInfo = DFUIsolatesInfo %>% select(Doern.lab.bank., patient_id)
DFUIsolatesInfo$DORN = paste0("DORN", DFUIsolatesInfo$Doern.lab.bank.)
DFUIsolatesInfo = DFUIsolatesInfo %>% left_join(CCmap, by="DORN")

CCmappingColors = CCmap %>% select(CCLabel) %>% unique() %>% left_join(ColorsCC,by="CCLabel")

UpdatedPhenotypes = read.csv("~/Desktop/GriceLabGit/Staphyloxanthin/Data/InVitroData/staphyloxanthin_paper_updated.csv")
UpdatedPhenotypes$DORN = sapply(UpdatedPhenotypes$IsolateID, function(x) str_replace(x, "SA", "DORN"))


# Reference genome the unitigs were mapped to
##############################################
USA300_LAC = read_tsv("data/dbgwas2022/FPR3757.gff", comment="#", col_names=F)
colnames(USA300_LAC) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr")

USA300_LAC %>% filter(feature=="region")
USA300_LAC = USA300_LAC %>% mutate(seqname = case_when(seqname=="contig1" ~ "NC_007793.1", 
                                                       seqname=="contig2" ~ "NC_007790.1", 
                                                       seqname=="contig3" ~ "NC_007791.1", 
                                                       seqname=="contig4" ~ "NC_007792.1"))

USA300_LAC = USA300_LAC %>% filter(feature=="CDS")

############################################################
# 1.  Making sense of top unitigs in KDE subset of Xanthin
############################################################
patternsSubsetXanthin = read_delim("data/dbgwas2022/XanthinSubsetZeroOne/step2/patterns.txt", " ",col_names=T)

UnitigsToPatterns = read.csv2("data/dbgwas2022/XanthinSubsetZeroOne/step1/gemma_input.unitig_to_pattern.binary",sep=" ", header=F)
colnames(UnitigsToPatterns) = c("unitigID", "pattern")
UnitigsToPatterns$NodeId = paste0("n",UnitigsToPatterns$unitigID )
UnitigsToPatterns = UnitigsToPatterns %>% select(NodeId, pattern)


# Staphyloxanthin (KDE clustered)
###############################
GWASResultsXanthinKDE =  data.frame(read_tsv("data/dbgwas2022/XanthinSubsetZeroOne/textualOutput/all_comps_nodes_info.tsv"))

bamfile = Rsamtools::scanBam("data/dbgwas2022/UnitigBams/StaphyloxanthinKDEUnitigs.bam")
NodeIDs = bamfile[[1]]$qname
positions = bamfile[[1]]$pos
contigs = bamfile[[1]]$rname
AlignmentFrameKDEXanthin = data.frame(NodeId=NodeIDs)
AlignmentFrameKDEXanthin$UnitigPosition= positions
AlignmentFrameKDEXanthin$ContigMapped = contigs
strands = bamfile[[1]]$strand
AlignmentFrameKDEXanthin$strand = strands
#AlignmentFrameKDEXanthin$Sequence =  bamfile[[1]]$seq

GWASResultsXanthinKDE = GWASResultsXanthinKDE %>% select(NodeId, AlleleFreq, q.Value, p.value, SequenceLength, Sequence) 

Sig10XanthinKDE = GWASResultsXanthinKDE %>% filter(q.Value < .1)

Sig10XanthinKDE = Sig10XanthinKDE %>% left_join(AlignmentFrameKDEXanthin, by="NodeId")
Sig10XanthinKDE$Annotation = NA
Sig10XanthinKDE$Overlap = NA
Sig10XanthinKDE$StartNext = NA
Sig10XanthinKDE$TotalGenes = NA

Annotations = c()
for(r in 1:nrow(Sig10XanthinKDE)){
  contigmapped = Sig10XanthinKDE[r, "ContigMapped"]
  if( !is.na(contigmapped) ){
    Sig10XanthinKDE[r, "Annotation"] = FindAnnotation(Sig10XanthinKDE[r, "UnitigPosition"],Sig10XanthinKDE[r, "SequenceLength"], USA300_LAC,Sig10XanthinKDE[r, "ContigMapped"],Sig10XanthinKDE[r,"strand"], 1)
    Sig10XanthinKDE[r, "Overlap"] = FindAnnotation(Sig10XanthinKDE[r, "UnitigPosition"],Sig10XanthinKDE[r, "SequenceLength"], USA300_LAC,Sig10XanthinKDE[r, "ContigMapped"], Sig10XanthinKDE[r,"strand"],2)
    Sig10XanthinKDE[r, "StartNext"] = FindAnnotation(Sig10XanthinKDE[r, "UnitigPosition"],Sig10XanthinKDE[r, "SequenceLength"], USA300_LAC,Sig10XanthinKDE[r, "ContigMapped"], Sig10XanthinKDE[r,"strand"],3)
    Sig10XanthinKDE[r, "TotalGenes"] = FindAnnotation(Sig10XanthinKDE[r, "UnitigPosition"],Sig10XanthinKDE[r, "SequenceLength"], USA300_LAC,Sig10XanthinKDE[r, "ContigMapped"], Sig10XanthinKDE[r,"strand"], 4)
    
  }
}
KDEpatternDist = read.csv2("data/dbgwas2022/XanthinSubsetZeroOne/step1/bugwas_input.all_rows.binary",sep=" ", header=T)



bamfile = NULL


# q values under .1:
#0.01161903 0.03842774 0.05393247 0.05470390 0.05643329 0.05947170 0.06219065 0.07189534 0.07192898 0.07954333


################
# q=0.01161903
###############
# 4 patterns: 722, 3979, 10, 663
patternsSubsetXanthin %>% filter(`q-value` <.02)

# Pattern 663:
###############
# pattern present in 4 genomes from 2 different patients
# Contains 168 different unitigs 
Nodes663 = UnitigsToPatterns %>% filter(pattern=="663")
Pattern663 = KDEpatternDist %>% filter(ps=="28419")
DFUIsolatesInfo %>% filter(DORN %in% colnames(Pattern663)[which(Pattern663==1)])


# Pattern 10:
###############
# pattern present in 4 genomes (3 from CC15, 4 from the same CC1-ST188 DFU)
Nodes10 =  UnitigsToPatterns %>% filter(pattern=="10")
Pattern10 = KDEpatternDist %>% filter(ps=="275355")
DFUIsolatesInfo %>% filter(DORN %in% colnames(Pattern10)[which(Pattern10==1)])
Sig10XanthinKDE %>% filter(NodeId %in% Nodes10$NodeId)
# n117232-- synonymous C to T in SA1027; non-synonymous (Y-->H) in SA1081; no transposon mutant
# n287325 -- synonymous G to A in SA1027; gene missing from SA1081
# n284323 -- upstream from yagU and downstream from argB
# n191066 -- in asnS; synonymous in 1027
# n267408 -- not in annotated region

# n275355 -- NS snp in SA1027; Nif3-like dinuclear metal center hexameric protein ; pgaptmp_001625

PresenceAbsence10 = data.frame(t(Pattern10[2:ncol(Pattern10)]))
colnames(PresenceAbsence10) = "AbsenceUnitig275355"
PresenceAbsence10$DORN = row.names(PresenceAbsence10)
PresenceAbsence10 = PresenceAbsence10 %>% left_join(CCmap,by="DORN")

PresenceAbsence10 = PresenceAbsence10 %>% left_join(UpdatedPhenotypes %>% select(DORN, staphyloxanthin))
colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence10$CCLabel))
PresenceAbsence_nif3_STX = ggplot(PresenceAbsence10, aes(x=factor(AbsenceUnitig275355), y=staphyloxanthin)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig275355), y=staphyloxanthin, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence_nif3_STX, file="data/dbgwas2022/UnitigPlots/n275355_NE1507.pdf", width=5,height=8)

# Pattern 3979:
###############
Nodes3979 =  UnitigsToPatterns %>% filter(pattern=="3979")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes3979$NodeId)
# One unitig (n56547) , located in uncharacterized protein which doesn't seem to be annotated by nebraska library

# Pattern 722:
###############
Nodes722 =  UnitigsToPatterns %>% filter(pattern=="722")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes722$NodeId)
# One unitig -- n108760 -- ebh 

pattern722 = KDEpatternDist %>% filter(ps=="108760")
PresenceAbsence722 = data.frame(t(pattern722[2:ncol(pattern722)]))
colnames(PresenceAbsence722) = "AbsenceUnitig108760"
PresenceAbsence722$DORN = row.names(PresenceAbsence722)
PresenceAbsence722 = PresenceAbsence722 %>% left_join(CCmap,by="DORN")

PresenceAbsence722 = PresenceAbsence722 %>% left_join(UpdatedPhenotypes %>% select(DORN, staphyloxanthin))
colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence722$CCLabel))
PresenceAbsence_Ebh_STX = ggplot(PresenceAbsence722, aes(x=factor(AbsenceUnitig108760), y=staphyloxanthin)) + geom_boxplot()+ geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig108760), y=staphyloxanthin, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence_Ebh_STX, file="data/dbgwas2022/UnitigPlots/n108760_NE1.pdf", width=5,height=8)


################
# q=0.03842774
###############
# 2 patterns: 4186, 1903
patternsSubsetXanthin %>% filter(`q-value` <.04 & `q-value`>.02)

# Pattern 4186
###############
Nodes4186 =  UnitigsToPatterns %>% filter(pattern=="4186")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes4186$NodeId)
Pattern4186 = KDEpatternDist %>% filter(ps=="223445")
DFUIsolatesInfo %>% filter(DORN %in% colnames(Pattern4186)[which(Pattern4186==1)])
# unitig n223445 -- ssbC ; missing from 6/110 isolates in 3 CCs
      # CC5s:  DORN333 (pgaptmp_000750)  & DORN339(pgaptmp_000202) -- synonymous
      # CC97: 1395 (pgaptmp_001273) -- synonymous
      # CC1: DORN1081 (pgaptmp_000675) , DORN1086, DORN1085 -- synonymous

# Pattern 1903
###############
Nodes1903 =  UnitigsToPatterns %>% filter(pattern=="1903")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes1903$NodeId)
Pattern1903 = KDEpatternDist %>% filter(ps=="72283")
DFUIsolatesInfo %>% filter(DORN %in% colnames(Pattern1903)[which(Pattern1903==1)])

    # DORN1086, DORN1085, DORN1081(pgaptmp_002164) -- synonymous
    # DORN1027 (pgaptmp_001685) -- NS (AA263 Thr->Ala)
    # DORN968 (missing)
    # DORN684 (pgaptmp_002223) -- NS (AA266 R->C)

PresenceAbsence1903  = data.frame(t(Pattern1903[2:ncol(Pattern1903)]))
colnames(PresenceAbsence1903) = "AbsenceUnitig72283"
PresenceAbsence1903$DORN = row.names(PresenceAbsence1903)
PresenceAbsence1903 = PresenceAbsence1903 %>% left_join(CCmap,by="DORN")

PresenceAbsence1903 = PresenceAbsence1903 %>% left_join(UpdatedPhenotypes %>% select(DORN, staphyloxanthin))
colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence1903$CCLabel))
PresenceAbsence_1903_STX = ggplot(PresenceAbsence1903, aes(x=factor(AbsenceUnitig72283), y=staphyloxanthin)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig72283), y=staphyloxanthin, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence_1903_STX, file="data/dbgwas2022/UnitigPlots/n72283_NE128.pdf", width=5,height=8)



################
# q=0.05393247
###############
patternsSubsetXanthin %>% filter(`q-value` <.054 & `q-value`  > .05)
# Pattern 661 only 
Nodes661 =  UnitigsToPatterns %>% filter(pattern=="661")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes661$NodeId)
# Both n266205 and n174975 map to proteins annotated as TIGR01741 family proteins(pgaptmp_000313 and pgaptmp_000310)
# present only in SA1086, 1081, 1395
# not possible to compare directly because of many closely related proteins in the genomes containing the unitig 
KDEpatternDist %>% filter(ps=="266205")


###############
# q=0.05470390
###############
patternsSubsetXanthin %>% filter(`q-value` <.055 & `q-value`  > .054)

# Pattern 311
##############
Nodes311 =  UnitigsToPatterns %>% filter(pattern=="311")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes311$NodeId)
# of the two unitigs, only n46811 is actually overlapping with a gene (pgaptmp_002617)
# Unitig is n46811 absent from :
#   DORN1027 (CC15) -- pgaptmp_001730 -- 100% identity
#   DORN186, DORN1085, DORN1081 (CC1-ST188) DORN1081 -- dont have a gene with more than 81% identity to this one
#   DORN1644 (CC5) -- Missing the gene

# Pattern 694
##############
Nodes694 =  UnitigsToPatterns %>% filter(pattern=="694")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes694$NodeId)
Pattern694 = KDEpatternDist %>% filter(ps=="30383")

PresenceAbsence694= data.frame(t(Pattern694[2:ncol(Pattern694)]))
colnames(PresenceAbsence694) = "AbsenceUnitigsPattern694"
PresenceAbsence694$DORN = row.names(PresenceAbsence694)
PresenceAbsence694 = PresenceAbsence694 %>% left_join(CCmap,by="DORN")

PresenceAbsence694 = PresenceAbsence694 %>% left_join(UpdatedPhenotypes %>% select(DORN, staphyloxanthin))
colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence694$CCLabel))
PresenceAbsence_694_STX = ggplot(PresenceAbsence694, aes(x=factor(AbsenceUnitigsPattern694), y=staphyloxanthin)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitigsPattern694), y=staphyloxanthin, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence_694_STX, file="data/dbgwas2022/UnitigPlots/n30383_NE1252.pdf", width=5,height=8)


# 4 unitigs missing from DORN1081, DORN1085, DORN1086, DORN1395, DORN1946
# n30383 -- phnE
# n50164 -- phnE
# n171021 -- phnE
# n284285 -- phnC
# PhnE:
  # DORN1081 -- pgaptmp_001856
      # n50164 --synonymous
      # n171021 -- 
      # n284285 -- 
  # DORN1395 -- pgaptmp_000098 
      # n50164 --synonymous
      # n171021 -- synonymous
      # n284285 -- synonymous
  # DORN1946 -- missing
# PhnC:
      # DORN1081(pgaptmp_001855) -- synonymous
      # DORN1395 (pgaptmp_000099) -- synonymous
      # DORN1946 (missing)
# NE952



###############
# q=0.05947170
###############
patternsSubsetXanthin %>% filter(`q-value` >.059 & `q-value`  <.06)

# Pattern 8365
##############
Nodes8365 =  UnitigsToPatterns %>% filter(pattern=="8365")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes8365$NodeId)
# aroK/shikimate kinase -- no transposon mutant

# Pattern 1677
##############
# Missing from 1027, 1081, 1085, 1086 and 807
# 1081: 
Nodes1677 =  UnitigsToPatterns %>% filter(pattern=="1677")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes1677$NodeId)
# n104738: gntK /gluconokinase
# 2637745	2639298
# gene missing from DORN807
# Many AA differences in SA1081, SA1085, SA1086
# Synonymous in SA1027

# Pattern 7397
##############
Nodes7397 =  UnitigsToPatterns %>% filter(pattern=="7397")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes7397$NodeId)
# Only one unitig (n301311), unmapped


# 0.07189534 0.07192898 0.07954333

###############
# q=0.06219065
###############
patternsSubsetXanthin %>% filter(`q-value` >.06 & `q-value`  <.07)

# Pattern 738
##############
Nodes738=  UnitigsToPatterns %>% filter(pattern=="738")
# Absent from  DORN1395, DORN1531, DORN1081, DORN1085,DORN1086 
Sig10XanthinKDE %>% filter(NodeId %in% Nodes738$NodeId)
# Only second two unitigs (n241250 and n225351) map to coding regions;
# n241250/pgaptmp_000325 (not annotated in transposon library)
# n225351/pgaptmp_000326 (pfoR in transposon library) /NE1254
# TAATAATACTGCAATAAGACGAGCTAATGGCGCTAAGATGACAATCGATCCAATTAAGTCG
# DORN1395 (pgaptmp_000272)-- Synonymous
# DORN1531 (nothing with significant identity)
# DORN1081(pgaptmp_001685) -- Synonymous


# Pattern 1869
###############
Nodes1869=  UnitigsToPatterns %>% filter(pattern=="1869")
# n244949
Sig10XanthinKDE %>% filter(NodeId %in% Nodes1869$NodeId)

Pattern1869 = KDEpatternDist %>% filter(ps=="244949")


# Unitig AATTTCCCAAGTATGGCACCTAAACCGAATAT
# pgaptmp_002647
# absent from 12 different genomes with low STX
  # DORN900 -- pgaptmp_002417 ; synonymous
  # DORN1395 -- pgaptmp_002435; synonymous
  # DORN1679 -- pgaptmp_001522 ; synonymous
  # DORN1974 -- pgaptmp_000594; synonymous
  # DORN915 -- pgaptmp_001903; synonymous
  # DORN962 -- pgaptmp_001717; synonymous
  # DORN968 -- pgaptmp_001478; synonymous
  # DORN1027 -- pgaptmp_001700; synonymous
  # DORN1081 -- pgaptmp_002179; synonymous
  # DORN1085 synonymous
  # DORN1086 synonymous
  # DORN807 -- missing gene
# Not a strong candidate given synonymous mutation, but widespread across CCs and transposon available
PresenceAbsence1869= data.frame(t(Pattern1869[2:ncol(Pattern1869)]))
colnames(PresenceAbsence1869) = "AbsenceUnitig244949"
PresenceAbsence1869$DORN = row.names(PresenceAbsence1869)
PresenceAbsence1869 = PresenceAbsence1869 %>% left_join(CCmap,by="DORN")

PresenceAbsence1869 = PresenceAbsence1869 %>% left_join(UpdatedPhenotypes %>% select(DORN, staphyloxanthin))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence1869$CCLabel))
PresenceAbsence1869_STX = ggplot(PresenceAbsence1869, aes(x=factor(AbsenceUnitig244949), y=staphyloxanthin)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig244949), y=staphyloxanthin, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence1869_STX, file="data/dbgwas2022/UnitigPlots/n244949_NE952.pdf", width=5,height=8)



###############
# q=0.07189534
###############

patternsSubsetXanthin %>% filter(`q-value` >0.071 & `q-value`  <.072)

# Pattern 1899
##############
Nodes1899=  UnitigsToPatterns %>% filter(pattern=="1899")
# 7 unitigs:
# n176102
# n355
# n90728
# n245935
# n240249
# n121979
# n301375
# Prioritized the 4 which are *present* in the minor allele genomes(DORN1086, DORN1085, DORN684, DORN1081): 
  # n355, n90728, n245935, n240249
  # 3 of these unmapped, one mapped nearby but not overlapping with coding sequence of a protein
  # blastn maps n90728 to  2277360-2277406 of DORN684,  2241080-2241126  of DORN1081 (fnbB in both these)
Pattern1899 = KDEpatternDist %>% filter(ps=="90728")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes1899$NodeId)

# Compare to DORN243 for CC45: Non-synonymous AA448 Q-->R comparead to 243
# Compare to DORN1000, DORN152, DORN1696 for CC1
# though there's a clear single AA difference in 243, this unitig corresponds to a whole variable stretch
# of the FnbB protein in other lineages (AA440-454 of DOR684's, with many variants compared to, say, DORN1000)


# Pattern 1780
###############
Nodes1780=  UnitigsToPatterns %>% filter(pattern=="1780")
# 9 different nodes, all 'absent' from the 5 MAF genomes
Pattern1899 = KDEpatternDist %>% filter(ps=="90728")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes1780$NodeId)

# Pattern 9
###########
Nodes9=  UnitigsToPatterns %>% filter(pattern=="9")
Sig10XanthinKDE %>% filter(NodeId %in% Nodes9$NodeId)
# both nodes unmapped and low allele frequency (3 genomes total; SA1081, SA1085, SA1027)
Pattern9= KDEpatternDist %>% filter(ps=="44396")


###############
# q=0.07954333
###############
patternsSubsetXanthin %>% filter(`q-value` >0.079 & `q-value`  <0.08)

# Pattern 8607
##############
Nodes8607=  UnitigsToPatterns %>% filter(pattern=="8607")
# 31 different unitigs 
Pattern8607 = KDEpatternDist %>% filter(ps=="203631")

# 12 of them map to a gene, and 5 of those are at allele frequency 6 (unitig *presence* is the MAF variant)
Subset8607 = Sig10XanthinKDE %>% filter(NodeId %in% Nodes8607$NodeId) %>% filter(!is.na(Overlap) & Overlap!=0) %>% filter(AlleleFreq==6)
# n48071-- pgaptmp_001539; hypothetical protein -- no transposon
# n199462 -- pgaptmp_001473; 1545736	1545921	
# n98200  pgaptmp_000318 – DUF4467 family protein which there are two of in FPR3757
# n134843 thrS (no transposon mutant)
# n203631 pgaptmp_000503 (532549..532929	mutant NE404 is right in the start of that)

PresenceAbsence8607= data.frame(t(Pattern8607[2:ncol(Pattern8607)]))
colnames(PresenceAbsence8607) = "PresenceUnitig203631"
PresenceAbsence8607$DORN = row.names(PresenceAbsence8607)
PresenceAbsence8607 = PresenceAbsence8607 %>% left_join(CCmap,by="DORN")

PresenceAbsence8607 = PresenceAbsence8607 %>% left_join(UpdatedPhenotypes %>% select(DORN, staphyloxanthin))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence8607$CCLabel))
PresenceAbsence8607_STX = ggplot(PresenceAbsence8607, aes(x=factor(PresenceUnitig203631), y=staphyloxanthin)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(PresenceUnitig203631), y=staphyloxanthin, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence8607_STX, file="data/dbgwas2022/UnitigPlots/n203631_NE404.pdf", width=5,height=8)


############################################
# 2. Making sense of the top biofilm unitigs
###########################################
# KDE-Subset biofilm
#####################

patternsSubsetBiofilm= read_delim("data/dbgwas2022/BiofilmSubsetZeroOne/step2/patterns.txt", " ",col_names=T)

GWASResultsBiofilmKDE =  data.frame(read_tsv("data/dbgwas2022/BiofilmSubsetZeroOne/textualOutput/all_comps_nodes_info.tsv"))
bamfile = Rsamtools::scanBam("data/dbgwas2022/UnitigBams/BiofilmKDEUnitigs.bam")
NodeIDs = bamfile[[1]]$qname
positions = bamfile[[1]]$pos
contigs = bamfile[[1]]$rname
strands=bamfile[[1]]$strand
AlignmentFrameKDEBiofilm = data.frame(NodeId=NodeIDs)
AlignmentFrameKDEBiofilm$UnitigPosition= positions
AlignmentFrameKDEBiofilm$ContigMapped = contigs
AlignmentFrameKDEBiofilm$strand=strands
#AlignmentFrameKDEBiofilm$Sequence =  bamfile[[1]]$seq


GWASResultsBiofilmKDE = GWASResultsBiofilmKDE %>% select(NodeId, AlleleFreq, q.Value, p.value, SequenceLength, Sequence) 

Sig10BiofilmKDE = GWASResultsBiofilmKDE %>% filter(q.Value < .1)

Sig10BiofilmKDE = Sig10BiofilmKDE %>% left_join(AlignmentFrameKDEBiofilm, by="NodeId")
Sig10BiofilmKDE$Annotation = NA
Sig10BiofilmKDE$Overlap = NA
Sig10BiofilmKDE$StartNext = NA
Sig10BiofilmKDE$TotalGenes = NA

Annotations = c()
for(r in 1:nrow(Sig10BiofilmKDE)){
  contigmapped = Sig10BiofilmKDE[r, "ContigMapped"]
  if( !is.na(contigmapped) ){
    Sig10BiofilmKDE[r, "Annotation"] = FindAnnotation(Sig10BiofilmKDE[r, "UnitigPosition"],Sig10BiofilmKDE[r, "SequenceLength"], USA300_LAC,Sig10BiofilmKDE[r, "ContigMapped"], Sig10BiofilmKDE[r, "strand"], 1)
    Sig10BiofilmKDE[r, "Overlap"] = FindAnnotation(Sig10BiofilmKDE[r, "UnitigPosition"],Sig10BiofilmKDE[r, "SequenceLength"], USA300_LAC,Sig10BiofilmKDE[r, "ContigMapped"],Sig10BiofilmKDE[r, "strand"], 2)
    Sig10BiofilmKDE[r, "StartNext"] = FindAnnotation(Sig10BiofilmKDE[r, "UnitigPosition"],Sig10BiofilmKDE[r, "SequenceLength"], USA300_LAC,Sig10BiofilmKDE[r, "ContigMapped"],Sig10BiofilmKDE[r, "strand"], 3)
    Sig10BiofilmKDE[r, "TotalGenes"] = FindAnnotation(Sig10BiofilmKDE[r, "UnitigPosition"],Sig10BiofilmKDE[r, "SequenceLength"], USA300_LAC,Sig10BiofilmKDE[r, "ContigMapped"], Sig10BiofilmKDE[r, "strand"], 4)
    
  }
}
write.csv(Sig10BiofilmKDE, file="data/dbgwas2022/BiofilmSubsetZeroOne/Sig10BiofilmKDE.csv")

UnitigsToPatterns = read.csv2("data/dbgwas2022/BiofilmSubsetZeroOne/step1/gemma_input.unitig_to_pattern.binary",sep=" ", header=F)
colnames(UnitigsToPatterns) = c("unitigID", "pattern")
UnitigsToPatterns$NodeId = paste0("n",UnitigsToPatterns$unitigID )
UnitigsToPatterns = UnitigsToPatterns %>% select(NodeId, pattern)



AlignmentFrameKDEBiofilm = AlignmentFrameKDEBiofilm %>% left_join(UnitigsToPatterns, by="NodeId")
AlignmentFrameKDEBiofilm = AlignmentFrameKDEBiofilm %>% left_join((patternsSubsetBiofilm %>% select(pattern, `q-value`, `p-value`)), by="pattern")
AlignmentFrameKDEBiofilm$UnitigPosition[is.na(AlignmentFrameKDEBiofilm$UnitigPosition)] <- 3000000
AlignmentFrameKDEBiofilm = AlignmentFrameKDEBiofilm %>% mutate(ContigMapped = if_else(is.na(ContigMapped), "Unmapped", as.character(ContigMapped)))
AlignmentFrameKDEBiofilm$NegLogTransformed = -log10(AlignmentFrameKDEBiofilm$`q-value`)
BiofilmKDENoPlasmids = AlignmentFrameKDEBiofilm %>% filter(ContigMapped %in% c("NC_007793.1", "Unmapped"))
KDEbiofilmPatternDist = read.csv2("data/dbgwas2022/BiofilmSubsetZeroOne/step1/bugwas_input.all_rows.binary", sep=" ")



# 8 different q values under .1 for biofilm
# 0.02958491 0.03328928 0.03491193 0.04624817 0.05085501 0.05368305 0.06323392 0.06449570
sort(unique(Sig10BiofilmKDE$q.Value))


###############
# q=0.02958491
###############


patternsSubsetBiofilm %>% filter(`q-value` < 0.03)
# 8 patterns 

# Pattern 16928
################
Nodes16928=  UnitigsToPatterns %>% filter(pattern=="16928")
# n94865
Sig10BiofilmKDE %>% filter(NodeId=="n94865")
# Absent from 34 different isolates from CC30, CC5, CC59, CC20, CC9, CC7
      # Absent from 14/33 CC5s, all CC7s, CC59s, CC9s
# Maps to AA482..AA500 pgaptmp_002780 (ABC transporter permease) 2766209-2768209, interrupted by NE1105 transposon 
Pattern16928= KDEbiofilmPatternDist %>% filter(ps=="94865")
DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern16928)[Pattern16928==1]))
# DORN781 (CC5 missing the unitig)-- pgaptmp_000127
# pgaptmp_000127 : D --> N in AA491


# Pattern 16913
################
Nodes16913=  UnitigsToPatterns %>% filter(pattern=="16913")
# Unmapped
Sig10BiofilmKDE %>% filter(NodeId=="n292821")


# Pattern 18125
################
Node18125 = UnitigsToPatterns %>% filter(pattern=="18125")
Pattern18125= KDEbiofilmPatternDist %>% filter(ps=="140188")
# missing from all 3 CC1-ST188 (Patient 310) isolates, 
# all CC20 isolates, all CC30 isolates, all CC45 isolates, all CC59 isolates, CC7, all CC9s and CC97; missing from a *subset* of CC5 isolates
# Present in all but 3 CC1s (the st188s), CC12, CC15, CC72, CC8

# n140188
Sig10BiofilmKDE %>% filter(NodeId=="n140188")
# Missing from 41 isolates from CC1,CC30,CC5,CC59,CC45,CC97,CC20,CC9,CC7
DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern18125)[Pattern18125==1]))
DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern18125)[Pattern18125==1]))
# pgaptmp_002477; 2460225-2461124	(interrupted by NE752)
# compare this gene 
# the unitig is present in DORN1373 (a CC5) but absent from DORN459
# DORN459-- pgaptmp_000490
# DORN1373 -- pgaptmp_002311
# But these two genes (DORN459 and 1373 version) are identical 


# Pattern 18074
################
Node18074 = UnitigsToPatterns %>% filter(pattern=="18074")
Pattern18074 = KDEbiofilmPatternDist %>% filter(ps=="215756")
Sig10BiofilmKDE %>% filter(NodeId=="n215756")
# Minor frequency pattern is the *absence* of unitig n215756
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern18074)[Pattern18074==1])) %>% select(CCLabel))
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern18074)[Pattern18074==1])) %>% select(CCLabel))
PresenceAbsence18074= data.frame(t(Pattern18074[2:ncol(Pattern18074)]))
colnames(PresenceAbsence18074) = "AbsenceUnitig215756"
PresenceAbsence18074$DORN = row.names(PresenceAbsence18074)
PresenceAbsence18074 = PresenceAbsence18074 %>% left_join(CCmap,by="DORN")
# This unitig, too, is present in a subset of CC5s, all CC72s, a subset of CC1s (the ST188's specifically), 
# all CC15s, CC7, and all CC8s
# But absent from CC12, a subset of CC1s, all CC20s, CC45s, CC59s, CC9s, and CC97, and a subset of CC5s
# Once again, not totally clonally correlated because of the CC5s
# pgaptmp_002772 (NE1068)
# unitig is absent from  DORN459-- pgaptmp_000779 (corresponds to M --> I in AA272 of DORN459)
# present in DORN1373 -- pgaptmp_002019

PresenceAbsence18074 = PresenceAbsence18074 %>% left_join(UpdatedPhenotypes %>% select(DORN, biofilm))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence18074$CCLabel))
PresenceAbsence18074_biofilm = ggplot(PresenceAbsence18074, aes(x=factor(AbsenceUnitig215756), y=biofilm)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig215756), y=biofilm, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence18074_biofilm, file="data/dbgwas2022/UnitigPlots/biofilm_n215756_NE1068.pdf", width=5,height=8)


# Pattern 4589
Node4589 = UnitigsToPatterns %>% filter(pattern=="4589")
Pattern4589 = KDEbiofilmPatternDist %>% filter(ps=="133456")
Sig10BiofilmKDE %>% filter(NodeId=="n133456")
# pgaptmp_000675 (dhaL; no transposon mutant)

# Pattern 5261
Node5261 = UnitigsToPatterns %>% filter(pattern=="5261")
Sig10BiofilmKDE %>% filter(NodeId=="n308928")
Pattern5261 = KDEbiofilmPatternDist %>% filter(ps=="308928")
# DUF1433 domain-containing protein (pgaptmp_002620); there are several 
# homologs to this gene in USA300 FPR3757 (with ~70% identity)
# corresponds to AA28-AA37 of pgaptmp_002620 (NE53 interrupts this)
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern5261)[Pattern5261==1])) %>% select(CCLabel))
(table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern5261)[Pattern5261==1])) %>% select(CCLabel)))
# again, a subset of CC5s is the only non-clonal grouping to be missing this unitig
# other than that: all CC1s, all CC12s, all CC15s, all CC45's, All CC7s, all CC72s, all CC8s, and all CC9s have it
# While all CC20s, all CC30s, 3 CC5s (DORN460, DORN468, DORN781), all CC59s, and all CC97s missing it
# comparing DORN1471's (pgaptmp_002212) to the unitig -- identical (synonymous)
# DORN767's  (pgaptmp_000376) -- S-->A at AA37 "DUF1433 domain-containing protein"
# DORN781's, DORN468's + 460's (all CC5s) missing the gene 
# DORN620 pgaptmp_002310
PresenceAbsence5621= data.frame(t(Pattern5261[2:ncol(Pattern5261)]))
colnames(PresenceAbsence5621) = "AbsenceUnitig308928"
PresenceAbsence5621$DORN = row.names(PresenceAbsence5621)
PresenceAbsence5621 = PresenceAbsence5621 %>% left_join(CCmap,by="DORN")

PresenceAbsence5621 = PresenceAbsence5621 %>% left_join(UpdatedPhenotypes %>% select(DORN, biofilm))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence5621$CCLabel))
PresenceAbsence5621_biofilm = ggplot(PresenceAbsence5621, aes(x=factor(AbsenceUnitig308928), y=biofilm)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig308928), y=biofilm, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence5621_biofilm, file="data/dbgwas2022/UnitigPlots/biofilm_n308928_NE53.pdf", width=5,height=8)





# Pattern 5264
#################
Node5264 = UnitigsToPatterns %>% filter(pattern=="5264")
Sig10BiofilmKDE %>% filter(NodeId=="n44463")
# 2608288	2608569 # tiny hypothetical protein interrupted by NE1175 (exclude bc it's not actually *in* this gene just between)
Pattern5264 = KDEbiofilmPatternDist %>% filter(ps=="44463")
# found upstream of the SSBP
# exculde 

# Minor allele frequency is hte absence of the unitig (23 isolates)
# CC20 CC30  CC5 CC59 CC97 
#   3    8    5    6    1 
# again, present/absent in a subset of CC5s but otherwise totally correlated with CC
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern5264)[Pattern5264==1])) %>% select(CCLabel))
DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern5264)[Pattern5264==1])) %>% arrange(CCLabel)
DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern5264)[Pattern5264==1])) %>% arrange(CCLabel)



# Pattern 5318
###############
Node5318 = UnitigsToPatterns %>% filter(pattern=="5318")
# Node n243947
Sig10BiofilmKDE %>% filter(NodeId=="n243947")
# Unitig *absent from*  24 genomes
  # including all CC20s, CC30s, CC97s; missing from 12 CC5s including:
      #DORN460 -- snp in ebh  -- G to D in AA738
      #DORN468 -- snp in ebh -- G to D in AA738
      #DORN781 -- snp in ebh -- G to D in AA738
    
# present in:DORN1881 608340...608372 (pgaptmp_000663)  though this gene is pseudogenized here (also in 1339)

Pattern5318 = KDEbiofilmPatternDist %>% filter(ps=="243947")

table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern5318)[Pattern5318==1])) %>% select(CCLabel))
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern5318)[Pattern5318==1])) %>% select(CCLabel))

PresenceAbsence5318= data.frame(t(Pattern5318[2:ncol(Pattern5318)]))
colnames(PresenceAbsence5318) = "AbsenceUnitig243947"
PresenceAbsence5318$DORN = row.names(PresenceAbsence5318)
PresenceAbsence5318 = PresenceAbsence5318 %>% left_join(CCmap,by="DORN")

PresenceAbsence5318 = PresenceAbsence5318 %>% left_join(UpdatedPhenotypes %>% select(DORN, biofilm))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence5318$CCLabel))
PresenceAbsence5318_biofilm = ggplot(PresenceAbsence5318, aes(x=factor(AbsenceUnitig243947), y=biofilm)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig243947), y=biofilm, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence5318_biofilm, file="data/dbgwas2022/UnitigPlots/biofilm_n243947_NE1.pdf", width=5,height=8)


###############
# q=0.03328928
###############
patternsSubsetBiofilm %>% filter(`q-value` < 0.034 & `q-value` > 0.033)

# Pattern 3845
Node3845 = UnitigsToPatterns %>% filter(pattern=="3845")
Sig10BiofilmKDE %>% filter(NodeId=="n289412")
# unmapped 

###############
# q=0.03491193
###############
patternsSubsetBiofilm %>% filter(`q-value` > 0.034 & `q-value` < 0.04)

# pattern 5263
##############
Node5263 = UnitigsToPatterns %>% filter(pattern=="5263")
Sig10BiofilmKDE %>% filter(NodeId=="n146611") 
Pattern5263 = KDEbiofilmPatternDist %>% filter(ps=="146611")

table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern5263)[Pattern5263==1])) %>% select(CCLabel))
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern5263)[Pattern5263==1])) %>% select(CCLabel))

PresenceAbsence5263= data.frame(t(Pattern5263[2:ncol(Pattern5263)]))
colnames(PresenceAbsence5263) = "AbsenceUnitig14611"
PresenceAbsence5263$DORN = row.names(PresenceAbsence5263)
PresenceAbsence5263 = PresenceAbsence5263 %>% left_join(CCmap,by="DORN")

PresenceAbsence5263 = PresenceAbsence5263 %>% left_join(UpdatedPhenotypes %>% select(DORN, biofilm))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence5263$CCLabel))
PresenceAbsence18074_biofilm = ggplot(PresenceAbsence5263, aes(x=factor(AbsenceUnitig14611), y=biofilm)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig14611), y=biofilm, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence18074_biofilm, file="data/dbgwas2022/UnitigPlots/biofilm_n14611_NE53.pdf", width=5,height=8)



setdiff(colnames(Pattern5263)[Pattern5263==1], colnames(Pattern5261)[Pattern5261==1])

# pattern 16914
###############
Node16914 = UnitigsToPatterns %>% filter(pattern=="16914")
Sig10BiofilmKDE %>% filter(NodeId=="n180922") # Doesn't map to an annotated gene

# Pattern 5266
##############
Node5266 = UnitigsToPatterns %>% filter(pattern=="5266")
Sig10BiofilmKDE %>% filter(NodeId=="n277928") # No transposon mutant for this


# Pattern 3385
##############

Node3385 = UnitigsToPatterns %>% filter(pattern=="3385")
Sig10BiofilmKDE %>% filter(NodeId=="n220775") # NE728
Pattern3385  = KDEbiofilmPatternDist %>% filter(ps=="220775")

# Unitig is absent from 17 genomes in 3 CCs
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern3385)[Pattern3385==1])) %>% select(CCLabel))
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern3385)[Pattern3385==1])) %>% select(CCLabel))
# DORN1455(CC5 missing unitig; pgaptmp_000032): A->V
# DORN1395 (CC97 missing unitig; pgaptmp_001911): synonymous 

PresenceAbsence3385= data.frame(t(Pattern3385[2:ncol(Pattern3385)]))
colnames(PresenceAbsence3385) = "AbsenceUnitig220775"
PresenceAbsence3385$DORN = row.names(PresenceAbsence3385)
PresenceAbsence3385 = PresenceAbsence3385 %>% left_join(CCmap,by="DORN")

PresenceAbsence3385 = PresenceAbsence3385 %>% left_join(UpdatedPhenotypes %>% select(DORN, biofilm))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence3385$CCLabel))
PresenceAbsence18074_biofilm = ggplot(PresenceAbsence3385, aes(x=factor(AbsenceUnitig220775), y=biofilm)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig220775), y=biofilm, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence18074_biofilm, file="data/dbgwas2022/UnitigPlots/biofilm_n220775_NE828.pdf", width=5,height=8)

###############
# q=0.04624817
###############
patternsSubsetBiofilm  %>% filter(`q-value` <.05 & `q-value`>.04)

# Pattern 907
##############
Node907 = UnitigsToPatterns %>% filter(pattern=="907")
Sig10BiofilmKDE %>% filter(NodeId=="n83471") 
# only missing from 2 genomes

# Pattern 1217
##############
Node1217= UnitigsToPatterns %>% filter(pattern=="1217")
Sig10BiofilmKDE %>% filter(NodeId=="n232380") 
Pattern1217  = KDEbiofilmPatternDist %>% filter(ps=="232380")

table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern1217)[Pattern1217==1])) %>% select(CCLabel))
table(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern1217)[Pattern1217==1])) %>% select(CCLabel))
#pgaptmp_000670 707016	707849 (NE1171 interrupts)
# Missing from one CC5 which is missing the gene (DORN2179)
# DORN1471 (CC20): pgaptmp_001342 -- synonymous
# DORN962 (CC59): pgaptmp_000809 -- synonymous
# DORN1395 (CC97): pgaptmp_000623 -- synonymous

# Pattern 5265
##############
Node5265= UnitigsToPatterns %>% filter(pattern=="5265")
Sig10BiofilmKDE %>% filter(NodeId=="n24226") 
# unmapped

# Pattern 5262
##############
Node5262= UnitigsToPatterns %>% filter(pattern=="5262")
Sig10BiofilmKDE %>% filter(NodeId=="n104087") 
Sig10BiofilmKDE %>% filter(NodeId=="n302888") 
# DUF1433 domain-containing protein again


#  0.05368305 0.06323392 0.06449570
###############
# q=0.05085501
###############
patternsSubsetBiofilm  %>% filter(`q-value` >.05 & `q-value`<.0509)

# Pattern 3389
##############
Node3389= UnitigsToPatterns %>% filter(pattern=="3389")
Sig10BiofilmKDE %>% filter(NodeId=="n231531") 
# unmapped


#   0.06323392 0.06449570
###############
# q=0.05368305
###############
patternsSubsetBiofilm  %>% filter(`q-value` >.053 & `q-value`<.054)

# Pattern 3843
##############
Node3843= UnitigsToPatterns %>% filter(pattern=="3843")
Sig10BiofilmKDE %>% filter(NodeId %in% Node3843$NodeId) 
# only 2 'variant' genomes for each of these unitigs

###############
# q=0.06323392
###############
patternsSubsetBiofilm  %>% filter(`q-value` >.06 & `q-value`<.064)

# Pattern 594
##############

Node594= UnitigsToPatterns %>% filter(pattern=="594")
Sig10BiofilmKDE %>% filter(NodeId %in% Node594$NodeId) 

Pattern594  = KDEbiofilmPatternDist %>% filter(ps=="5886")

PresenceAbsence594= data.frame(t(Pattern594[2:ncol(Pattern594)]))
colnames(PresenceAbsence594) = "AbsenceUnitig144640_177054"
PresenceAbsence594$DORN = row.names(PresenceAbsence594)
PresenceAbsence594 = PresenceAbsence594 %>% left_join(CCmap,by="DORN")

PresenceAbsence594 = PresenceAbsence594 %>% left_join(UpdatedPhenotypes %>% select(DORN, biofilm))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence594$CCLabel))
PresenceAbsence594_biofilm = ggplot(PresenceAbsence594, aes(x=factor(AbsenceUnitig144640_177054), y=biofilm)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig144640_177054), y=biofilm, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence594_biofilm, file="data/dbgwas2022/UnitigPlots/biofilm_n144640_n177054_NE216.pdf", width=5,height=8)


# 2 unitigs in dhaK and 1 in dhaL; missing from 5 (CC97, one CC5, one CC20)
(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern594)[Pattern594==1])) %>% arrange(CCLabel))
#dhaKs:
# GGTGGTACGCCGTTATCTGAATTAAATATCGTAACTAAATATATTCAACAAAATTT -- corresponds to AA255 to AA272 of FPR3757
# AATGTTGCTAAATGGTTTGTTGGTGATTATATGACATCTTTAGACATGCAAGGTT -- corresponds to AA279 to AA296 of FPR3757


# DORN1471 (CC20 missing unitig): pgaptmp_001338 ; n144640 corresponds to I->V AA265;n177054 corresponds to N-->D AA287
# DORN2179 (CC5 missing unitig) : missing gene
# DORN1395 (CC97 missing unitig): pgaptmp_000627 ; n144640 corresponds to I->V AA265; n177054 corresponds to N-->D AA287


###############
# q=0.06449570
###############
patternsSubsetBiofilm  %>% filter(`q-value` >.064 & `q-value`<.065)

# Pattern 3301
Node3301= UnitigsToPatterns %>% filter(pattern=="3301")
Sig10BiofilmKDE %>% filter(NodeId %in% Node3301$NodeId) #bshA 1516184	1517326 -- NE1728
Pattern3301  = KDEbiofilmPatternDist %>% filter(ps=="269297")
(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern3301)[Pattern3301==1])) %>% arrange(CCLabel))
(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern3301)[Pattern3301==1])) %>% arrange(CCLabel))
# Missing from CC20s, subset of CC5s, CC97s
# Missing from DORN1455(CC5)--pgaptmp_000627 -- F->Y AA68
# Missing from DORN1471 (CC20)--pgaptmp_000566 --synonymous
# Missing from DORN1395  (CC97)--pgaptmp_001386 --synonymous









# Pattern 6302
Node6302= UnitigsToPatterns %>% filter(pattern=="6302")
Sig10BiofilmKDE %>% filter(NodeId %in% Node6302$NodeId) #pgaptmp_002730; 2718911	2719615 ; interrupted by NE1527
Pattern6302  = KDEbiofilmPatternDist %>% filter(ps=="221501")
(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter((DORN %in% colnames(Pattern6302)[Pattern6302==1])) %>% arrange(CCLabel))
(DFUIsolatesInfo %>% filter(DORN %in% colnames(KDEbiofilmPatternDist) ) %>% filter(!(DORN %in% colnames(Pattern6302)[Pattern6302==1])) %>% arrange(CCLabel))

# spans AA85-AA99 of pgaptmp_002730 in FPR3757 
# absent from exactly one CC15/5 of them in the comparison
# absent from all CC30s, subset of 14 CC5s, all CC20s, all CC59s

# Absent from:
# DORN1027 (CC15) -- pgaptmp_001622 -- Y-> S at AA91
# DORN1471 (CC20) --pgaptmp_002103 -- T -> A at AA90
# DORN1410 (CC30) -- pgaptmp_002129 -- T->A at AA90
# DORN283 (CC5) -- pgaptmp_002500 -- P-->T at AA95
# DORN900 (CC59) -- pgaptmp_002504 -- Synonymous 
PresenceAbsence6302= data.frame(t(Pattern6302[2:ncol(Pattern6302)]))
colnames(PresenceAbsence6302) = "AbsenceUnitig221501"
PresenceAbsence6302$DORN = row.names(PresenceAbsence6302)
PresenceAbsence6302 = PresenceAbsence6302 %>% left_join(CCmap,by="DORN")

PresenceAbsence6302 = PresenceAbsence6302 %>% left_join(UpdatedPhenotypes %>% select(DORN, biofilm))

colors_To_use = CCmappingColors %>% arrange(CCLabel) %>% filter(CCLabel %in% unique(PresenceAbsence6302$CCLabel))
PresenceAbsence6302_biofilm = ggplot(PresenceAbsence6302, aes(x=factor(AbsenceUnitig221501), y=biofilm)) + geom_boxplot()+
  geom_jitter(height=0, width=.2, aes(x=factor(AbsenceUnitig221501), y=biofilm, color=CCLabel)) + scale_color_manual(values=colors_To_use$hexval) + theme_classic() + stat_compare_means()
ggsave(PresenceAbsence6302_biofilm, file="data/dbgwas2022/UnitigPlots/biofilm_n221501_NE1527.pdf", width=5,height=8)








