# Amy Campbell
# April 2022
# Taking in SAM alignment from BWA Mem of dbgwas-generated unitigs against FPRP3757

# Required packages
###################
library(tidyverse)
library(dplyr)
library(Rsamtools)


# Functions
###########

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

# Annotation file

FPR3727 = read_tsv("Documents/DFUData/GWAS/USA300_FPR3757_NoSeqs.gff", skip=5, col_names=F) #[1:2642, ]
colnames(FPR3727) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr")


# Read in DBGWAS Tabular Output
###############################
NodeInfo98 = read.csv2("Documents/DFUData/GWAS/all_comps_nodes_info98.tsv",sep='\t',row.names=as.character(1:7683))
NodeInfo104 = read.csv2("Documents/DFUData/GWAS/all_comps_nodes_info104.tsv",sep='\t',row.names=as.character(1:10618))
NodeInfo219 = read.csv2("Documents/DFUData/GWAS/all_comps_nodes_info219.tsv",sep='\t',row.names=as.character(1:3169))

NodeInfo98 = NodeInfo98[,2:17]
NodeInfo104 = NodeInfo104[,2:17]
NodeInfo219 = NodeInfo219[,2:17]
ColNamesVect = c("NodeID", "AlleleFreq", "Pheno0Count","Pheno0TotalCount","Pheno1Count", "Pheno1TotalCount","NACount", "NATotalCount", "Significant", "pvalue", "qValue", "EstEffect", "WaldStat", "Sequence", "SequenceLength", "Annotations")


colnames(NodeInfo98)= ColNamesVect
colnames(NodeInfo104)= ColNamesVect
colnames(NodeInfo219)= ColNamesVect

# Read in BamFiles
#################
bamfile98 = Rsamtools::scanBam("Documents/DFUData/GWAS/UnitigAlignment98.bam")
bamfile104 = Rsamtools::scanBam("Documents/DFUData/GWAS/UnitigAlignment104.bam")
bamfile219 = Rsamtools::scanBam("Documents/DFUData/GWAS/UnitigAlignment219.bam")

# Set up 98 cluster
NodeIDs98 = bamfile98[[1]]$qname
positions98= bamfile98[[1]]$pos
AlignmentFrame98 = data.frame(rep(0, length(NodeIDs98)))
AlignmentFrame98$NodeID = NodeIDs98  
AlignmentFrame98$UnitigPosition= positions98

AlignmentFrame98 = AlignmentFrame98 %>% left_join((NodeInfo98 %>% select(NodeID, Sequence, Annotations,SequenceLength, qValue, pvalue)), by="NodeID")



max(AlignmentFrame98$UnitigPosition, na.rm=T)
AlignmentFrame98$UnitigPosition <- AlignmentFrame98$UnitigPosition %>% replace_na(3000000)

# Set up 104 clusters
NodeIDs104 = bamfile104[[1]]$qname
positions104= bamfile104[[1]]$pos

AlignmentFrame104 = data.frame(rep(0, length(NodeIDs104)))
AlignmentFrame104$NodeID = NodeIDs104
AlignmentFrame104$UnitigPosition= positions104

AlignmentFrame104 = AlignmentFrame104 %>% left_join((NodeInfo104 %>% select(NodeID, Sequence, Annotations,SequenceLength, qValue, pvalue)), by="NodeID")

# Set up 219 clusters
NodeIDs219 = bamfile219[[1]]$qname
positions219= bamfile219[[1]]$pos

AlignmentFrame219 = data.frame(rep(0, length(NodeIDs219)))
AlignmentFrame219$NodeID = NodeIDs219
AlignmentFrame219$UnitigPosition= positions219


AlignmentFrame219 = AlignmentFrame219 %>% left_join((NodeInfo219 %>% select(NodeID, Sequence, Annotations,SequenceLength, qValue, pvalue)), by="NodeID")


max(AlignmentFrame219$UnitigPosition, na.rm=T)
AlignmentFrame219$UnitigPosition <- AlignmentFrame219$UnitigPosition %>% replace_na(3000000)
#AlignmentFrame219 = AlignmentFrame219 %>% select(NodeID, UnitigPosition)



FindAnnotation(FPR3727, 1)


FindAnnotation(677362,58, FPR3727, 1)  
FindAnnotation(677362,58, FPR3727, 2)  


StaphyloxanthinPaperData=read.csv("Documents/DFUData/GWAS/staphyloxanthin_paper_data.csv")
StaphyloxanthinPaperData$DORN = paste0("DORN", StaphyloxanthinPaperData$DORN)
StaphyloxanthinPaperData= StaphyloxanthinPaperData %>% select(DORN, staphyloxanthin)
for(r in 1:nrow(AlignmentFrame219)){
  row = AlignmentFrame219[r,]
  AlignmentFrame219[r, "Annotation3757"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 1)
  AlignmentFrame219[r, "AnnotationOverlap"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 2)
  AlignmentFrame219[r, "AnnotationLocus"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 3)
}

write.csv(AlignmentFrame219, file="Documents/DFUData/GWAS/AnnotatedUnitigs219.csv", quote=F)


for(r in 1:nrow(AlignmentFrame98)){
  row = AlignmentFrame98[r,]
  AlignmentFrame98[r, "Annotation3757"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 1)
  AlignmentFrame98[r, "AnnotationOverlap"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 2)
  AlignmentFrame98[r, "AnnotationLocus"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 3)
}

write.csv(AlignmentFrame98, file="Documents/DFUData/GWAS/AnnotatedUnitigs98.csv", quote=F)






max(AlignmentFrame104$UnitigPosition, na.rm=T)
AlignmentFrame104$UnitigPosition <- AlignmentFrame104$UnitigPosition %>% replace_na(3000000)

for(r in 1:nrow(AlignmentFrame104)){
  row = AlignmentFrame104[r,]
  AlignmentFrame104[r, "Annotation3757"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 1)
  AlignmentFrame104[r, "AnnotationOverlap"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 2)
  AlignmentFrame104[r, "AnnotationLocus"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], FPR3727, 3)
}
write.csv(AlignmentFrame104, file="Documents/DFUData/GWAS/AnnotatedUnitigs104.csv", quote=F)


#TTAAGTTATTAGTGCCTCTTATGCAGTTGCTCA
Unitigs_By_Isolate_Path219 = "Documents/DFUData/GWAS/bugwas_input.all_rows219.binary"
UnitigDist219 = read.table(Unitigs_By_Isolate_Path219,header = T )
UnitigDist219$NodeID = paste0("n", UnitigDist219$ps)

Node261087_distribution = UnitigDist219 %>% filter (NodeID =="n261087")
NodeInfo219 %>% filter(NodeID=="n261087")


Unitigs_By_Isolate_Path104 = "Documents/DFUData/GWAS/bugwas_input.all_rows104.binary"
UnitigDist104 = read.table(Unitigs_By_Isolate_Path104,header = T )
UnitigDist104$NodeID = paste0("n", UnitigDist104$ps)

# DORN900 has it 

# Go through all the unitigs shared between the 104 subset and the 219 full set
sharedsequencelist = c("TGAAAAAATACGCGGCACAAGGTAAGTTTGCGGAAG",
"TATTTTCACAAAATACTATAATGAGGATAGTAAATAGAGAGGAG",
"TTGAAAAAATAGAAGCCAACATTAGATAATTCAATGAAATATG",
"AATATTTAGTAAACATGGTGAACAATATTTCAGGAATT",
"AATTTCCCAAGTATGGCACCTAAACCGAATAT",
"ATTCGGTTTAGGTGCCATACTTGGGAAATTA", 
"TTAAGTTATTAGTGCCTCTTATGCAGTTGCTCA")

SharedannotationList=c("arcC", "YagU (promoter region)", "SSL5","Shikimate kinase", "gluconate:H+ symporter (1)","gluconate:H+ symporter (2)", "Unknown" )
plots219=c()
plots104=c()

i=0
for(sequencestring in sharedsequencelist){
  i=i+1
  nodeID_219=(AlignmentFrame219 %>% filter(Sequence==sequencestring))$NodeID
  DistributionNode219 = UnitigDist219 %>% filter(NodeID == nodeID_219)
  DistributionNode219$ps=NULL
  DistributionNode219$NodeID=NULL
  print(SharedannotationList[[i]])
  print(DistributionNode219)
  
  NodePresent = data.frame(t(DistributionNode219))
  colnames(NodePresent) = c("UnitigPresent")
  NodePresent$DORN = row.names(NodePresent)
  NodePresent =  NodePresent %>% left_join(StaphyloxanthinPaperData, by="DORN")
  print(sum(NodePresent$UnitigPresent))
  
  plots219[[i]] = ggplot(NodePresent, aes(x=factor(UnitigPresent), y=staphyloxanthin, group=factor(UnitigPresent))) + geom_boxplot(fill="goldenrod") + ggtitle(paste0("Staphyloxanthin Phenotype for 219 Isolates: ", SharedannotationList[[i]]))
  print(t.test(staphyloxanthin~UnitigPresent,data=NodePresent))
  
  
  nodeID_104=(AlignmentFrame104 %>% filter(Sequence==sequencestring))$NodeID
  DistributionNode104 = UnitigDist104 %>% filter(NodeID == nodeID_104)
  DistributionNode104$ps=NULL
  DistributionNode104$NodeID=NULL
  NodePresent = data.frame(t(DistributionNode104))
  colnames(NodePresent) = c("UnitigPresent")
  NodePresent$DORN = row.names(NodePresent)
  NodePresent =  NodePresent %>% left_join(StaphyloxanthinPaperData, by="DORN")
  print(sum(NodePresent$UnitigPresent))
  plots104[[i]] = ggplot(NodePresent, aes(x=factor(UnitigPresent), y=staphyloxanthin, group=factor(UnitigPresent))) + geom_boxplot(fill="goldenrod") + ggtitle(paste0("Staphyloxanthin Phenotype for 104 Isolates: ", SharedannotationList[[i]]))
  print(t.test(staphyloxanthin~UnitigPresent,data=NodePresent))
  
    }


ggsave(gridExtra::grid.arrange(plots219[[1]], plots104[[1]],ncol=2), file="arcC_Shared_Unitigs_Staphyloxanthin.png", height=6, width=13)
ggsave(gridExtra::grid.arrange(plots219[[2]], plots104[[2]],ncol=2), file="yagU_Shared_Unitigs_Staphyloxanthin.png", height=6, width=13)
ggsave(gridExtra::grid.arrange(plots219[[3]], plots104[[3]],ncol=2), file="ssl5_Shared_Unitigs_Staphyloxanthin.png", height=6, width=13)
ggsave(gridExtra::grid.arrange(plots219[[4]], plots104[[4]],ncol=2), file="shikimate_Shared_Unitigs_Staphyloxanthin.png", height=6, width=13)
ggsave(gridExtra::grid.arrange(plots219[[5]], plots104[[5]],ncol=2), file="gluconate_Shared_Unitigs_Staphyloxanthin.png", height=6, width=13)
ggsave(gridExtra::grid.arrange(plots219[[6]], plots104[[6]],ncol=2), file="gluconate2_Shared_Unitigs_Staphyloxanthin.png", height=6, width=13)
ggsave(gridExtra::grid.arrange(plots219[[7]], plots104[[7]],ncol=2), file="Unknown_Shared_Unitigs_Staphyloxanthin.png", height=6, width=13)


NodeInfo219 %>% filter(NodeID=="n8782")
# ArcC unitig:
# Present in DORN1471 and DORN1395 (in 219 isolate set) , absent from DORN1540

# YagU unitig:
# present in 14 (in 219 isolate set) including DORN962, absent from DORN1471
# Actually *absent* from dorn962, present in DORN1471


# SSL5 unitig 
# Present in 58 (in 219 isolate set)  including DORN999, absent from DORN1471
# actually *absent* from DORN999, present in DORN1471

# Shikimate kinase unitig:
# Absent from  45 (in 219 isolate set)  including DORN105, DORN657, DORN892,DORN1502
# Present in  DORN999

#gluconate:H+ symporter
# Absent from 16 (in 219 isolates set) including DORN962, DORN1081; present in DORN1340, DORN892, DORN1502

# gluconate:H+ symporter # 2 unitig:
# Absent from 54 including DORN962, DORN1502; present in DORN1197, DORN881

# unknown unitig:
# Present in 12 including 1081, 962; absent from DORN1471, DORN1353

PresencedataGluconateSymporter2=(plots219[[6]])$data
PresencedataGluconateSymporter1=(plots219[[5]])$data

PresenceGluconate1 = PresencedataGluconateSymporter1 %>% filter(UnitigPresent==1)

PresenceGluconate2 = PresencedataGluconateSymporter2 %>% filter(UnitigPresent==1)

# Every DORN containing gluconate sypmorter 1 unitig also present in gluconate symporter 2
setdiff(PresenceGluconate1$DORN, PresenceGluconate2$DORN)
setdiff(PresenceGluconate2$DORN, PresenceGluconate1$DORN)


