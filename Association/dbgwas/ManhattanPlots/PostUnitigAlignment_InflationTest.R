# Amy Campbell
# January 2022
# Taking in SAM alignment from BWA Mem of dbgwas-generated unitigs against SA_502A 
# Specifically for the inflation test (out of curiosity to see if the 1 or 2 hits were in the same part of the genome)

# Required packages
###################
library(tidyverse)
library(dplyr)
library(Rsamtools)
library(qqman)

# Read in DBGWAS Tabular Output
###############################
setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")
NodeInfo = data.frame(read_tsv("data/dbgwasInflationOutput/textualOutput/all_comps_nodes_info.tsv"))
NodeInfoFull = NodeInfo
NodeInfo$Annotations.sep..... = NULL
bamfile = Rsamtools::scanBam("data/UnitigsInflationTest.bam")

NodeIDs = bamfile[[1]]$qname
positions = bamfile[[1]]$pos


AlignmentFrame = data.frame(rep(0, length(NodeIDs)))
AlignmentFrame$NodeId = NodeIDs   
AlignmentFrame$UnitigPosition= positions


UnalignedUnitigs = AlignmentFrame %>% filter(is.na(UnitigPosition))
UnalignedUnitigsInfo = NodeInfo %>% filter(NodeId %in% UnalignedUnitigs$NodeId)

# SA502A genome 
###############
# Chromosome
# 1:2764699
# Plasmid 
# 2764700:2787689

# 257/ 488 of the unitigs that are significant at q-value <.01 level are not mapped :( 
# Change positions that are NA to 3000000 which i'll label as 'non-mapping' 
max(AlignmentFrame$UnitigPosition, na.rm=T)
AlignmentFrame$UnitigPosition <- AlignmentFrame$UnitigPosition %>% replace_na(3000000)
AlignmentFrame = AlignmentFrame %>% select(NodeId, UnitigPosition)

# Merge to attach unitig position to node info 
GWAS_result = NodeInfo %>% select(NodeId, AlleleFreq, q.Value, p.value, SequenceLength) %>% left_join(AlignmentFrame, by="NodeId")
hist(GWAS_result$UnitigPosition)
#GWAS_result = GWAS_result %>% filter(AlleleFreq < 110)
# GWAS_result = GWAS_result %>% mutate(Chromosome = case_when(UnitigPosition<=2764699 ~ "SA502A Chromosome", 
#                                               (UnitigPosition > 2764699 & UnitigPosition<3000000) ~ "Unnamed Plasmid", 
#                                               UnitigPosition==3000000 ~ "Unmapped"))
# 
# GWAS_result = GWAS_result %>% mutate(ChromosomeNumeric = case_when(UnitigPosition<3000000 ~ 1, 
#                                                             (UnitigPosition >=3000000 ) ~ 2)) 
#                                                 
GWAS_result = GWAS_result %>% filter(!is.na(q.Value))
GWAS_result$logtransformedQ = log10(GWAS_result$q.Value)
# GWAS_result$ChromosomeNumeric=1
# 
# 
# ggplot(GWAS_result, aes(x=UnitigPosition, y=-logtransformedQ)) + geom_point() + xlab("Unitig Position on S. aureus 502A") + ylab("-Log-10-transformed Q Value") + geom_line(y=-log10(.01))
# 
# qq(GWAS_result$q.Value)




unmapped = GWAS_result %>% filter(UnitigPosition == 3000000)



# Filter to variants that are in <208 and > 11 isolates (.95*219=208.05)
########################################################################
GWAS_resultFilteredAlleleFreq = GWAS_result 
ggplot(GWAS_result, aes(x=UnitigPosition, y=-logtransformedQ)) +
  geom_point() + xlab("Unitig Position on S. aureus 502A") +
  ylab("-Log-10-transformed Q Value") + geom_line(y=-log10(.01), linetype="dashed", colour="red") +
  geom_line(y=-log10(.05), linetype="dashed", colour="blue") #+ geom_text()

qq(GWAS_resultFilteredAlleleFreq$p.value)

# Looking at distributions of these variants
############################################
significant05= (GWAS_resultFilteredAlleleFreq %>% filter(q.Value < .05))
View(significant05)
View(NodeInfoFull %>% filter(NodeId %in% significant05$NodeId))

SA502A_annotations = (read_tsv("data/SA_502A.gff", skip=3, col_names=F))[1:2642, ]
colnames(SA502A_annotations) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr")

FindAnnotation <- function(unitigposition,unitiglength, AnnotationFrame, returnType){
  # INPUTS:
    # - unitigposition - starting position of unitig 
    # - unitiglength - sequence length of the unitig in question
    # - AnnotationFrame - the GFF data frame of the annotated reference genome 
    # - returnType
      # 1 : Annotation value; tells the annotation value
      # 2 : Overlaps ; tells how many bp this unitig overlaps with an actual CDS annotation by
      # 3:  Tells the start locus of the closest annotation (useful in cases of no overlap)
  # OUTPUTS:
    # Either full annotation string from .gff located at the CDS start closest to the unitig start
    # or the # of overlapping basepairs between the annotation and unitig depending on returnType val

  # Can't look for unitig that didn't map anywhere

  unitigposition = (as.integer(unitigposition))
  unitiglength = as.integer(unitiglength)
  print(unitigposition)
  if(unitigposition==3000000){
  return(NA)}
  
  # Otherwise, find row unitig is closest to 
  closest_start_Unitig_start = which.min( abs(AnnotationFrame$start-unitigposition) )
  Overlap = (intersect(  ((AnnotationFrame[closest_start_Unitig_start,])$start : (AnnotationFrame[closest_start_Unitig_start,])$end), unitigposition:(unitigposition + unitiglength)))
  if(returnType == 1){
    return((AnnotationFrame[closest_start_Unitig_start,])$attr)
  }else if(returnType == 2){
    return(length(Overlap))
  }else if(returnType==3){
    return((AnnotationFrame[closest_start_Unitig_start,])$start)
  }
  }

significant05 = significant05 %>% filter(!is.na(UnitigPosition))
for(r in 1:nrow(significant05)){
  row = significant05[r,]
  significant05[r, "Annotation502A"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 1)
  significant05[r, "AnnotationOverlap"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 2)
  significant05[r, "AnnotationLocus"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 3)
}

