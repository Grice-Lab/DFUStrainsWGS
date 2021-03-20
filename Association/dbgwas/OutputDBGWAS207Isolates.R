library(tidyr)
library(readr)
library(ggplot2)
library(dplyr)
library(ggtree)
library(ape)

setwd("~/Desktop/GriceLabGit/DFUStrainsWGS")

# Define file paths 
####################
patterns_bugwas = "/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning_207Isolates/DBGWAS_Output_207_3-17/step2/patterns.txt"
NodeInfo = data.frame(read_tsv("/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning_207Isolates/DBGWAS_Output_207_3-17/textualOutput/all_comps_nodes_info.tsv"))
Alignment502A = "/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning_207Isolates/DBGWAS_Output_207_3-17/UnitigAlignment207.bam"
UnitigsToPatterns  = "/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning_207Isolates/DBGWAS_Output_207_3-17/step1/gemma_input.unitig_to_pattern.binary"
Annotations502A = "data/SA_502A.gff"


# Hacking qqman  qq() fxn to plot qq plots 
#################################################################################################
bugWASpatterns = read_delim(patterns_bugwas, " ", col_names=T)

pvect =  bugWASpatterns$`p-value`
qqobserved = -log10(sort(pvect,decreasing=FALSE))
qqexpect = -log10( ppoints(length(pvect) ))

dataframe = data.frame(expect=qqexpect, observe=qqobserved)
ggplot(data=dataframe, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for DBGWAS p-values")  + theme_classic()#+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)

# Alignment of all Unitigs to SA502A
########################################
bamfile = Rsamtools::scanBam(Alignment502A)
NodeIDs = bamfile[[1]]$qname
positions = bamfile[[1]]$pos

AlignmentFrame = data.frame(rep(0, length(NodeIDs)))
AlignmentFrame$NodeId = NodeIDs   
AlignmentFrame$UnitigPosition= positions


# Unitigs to bugWAS patterns mapping
unitig_pattern = read_delim(UnitigsToPatterns, delim=" ", col_names=F)
colnames(unitig_pattern) = c("NodeId","pattern")
unitig_pattern$NodeId = paste0("n", unitig_pattern$NodeId)

# Filter unitig_pattern mapping to include only the patterns that were tested by bugWAS
unitig_pattern = unitig_pattern %>% filter(pattern %in% bugWASpatterns$pattern)

# also filter AlignmentFrame
AlignmentFrame = AlignmentFrame %>% filter(NodeId %in% unitig_pattern$NodeId)


# SA502A genome 
###############
# Chromosome
# 1:2764699
# Plasmid 
# 2764700:2787689
max(AlignmentFrame$UnitigPosition, na.rm=T)
# Maximum position anything aligned to was 2764624, so 3000000 is a safe 'didn't align' position 

AlignmentFrame$UnitigPosition <- AlignmentFrame$UnitigPosition %>% replace_na(3000000)

AlignmentFrame = AlignmentFrame %>% left_join(unitig_pattern, by="NodeId")

colnames(bugWASpatterns) = c("pattern", "pvalue", "qvalue", "weight", "wald_statistic")

GWAS_result = bugWASpatterns %>% select(pattern, pvalue, qvalue) %>% left_join(AlignmentFrame, by="pattern")

GWAS_result = GWAS_result %>% filter(!is.na(qvalue))

GWAS_result$logtransformedQ =  log10(GWAS_result$qvalue)

manhattanplot <- ggplot(GWAS_result, aes(x=UnitigPosition, y=-logtransformedQ)) +
  geom_point() + xlab("Unitig Position on S. aureus 502A") +
  ylab("-Log10(Q-Value)") + geom_line(y=-log10(.01), linetype="dashed", colour="red") +
  geom_line(y=-log10(.05), linetype="dashed", colour="blue")+ theme_classic() +theme(plot.title=element_text(size=20, hjust=.5)) + annotate("text", label="q-value=.01", color="red", x= 100000, y=2.1) + annotate("text", label="q-value=.05", color="blue", x= 100000, y=1.4) + ggtitle("DBGWAS Unitigs Mapped to SA502A")

ggsave(manhattanplot, file="manhattanplot.pdf", width=15, height=10)


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
  
  if(unitigposition==3000000){
    return(NA)}
  
  # Otherwise, find row unitig is closest to 
  closest_start_Unitig_start = which.min( abs(AnnotationFrame$start-unitigposition) )
  Overlap = (intersect( ((AnnotationFrame[closest_start_Unitig_start,])$start : (AnnotationFrame[closest_start_Unitig_start,])$end), unitigposition:(unitigposition + unitiglength)))
  if(returnType == 1){
    return((AnnotationFrame[closest_start_Unitig_start,])$attr)
  }else if(returnType == 2){
    return(length(Overlap))
  }else if(returnType==3){
    return((AnnotationFrame[closest_start_Unitig_start,])$start)
  }
}
significant01= (GWAS_result %>% filter(qvalue < .01))
significant01 = significant01 %>% left_join(NodeInfo, by = "NodeId")

for(r in 1:nrow(significant01)){
  row = significant01[r,]
  significant01[r, "Annotation502A"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 1)
  significant01[r, "AnnotationOverlap"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 2)
  significant01[r, "AnnotationLocus"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 3)
}
# Interpreting 
