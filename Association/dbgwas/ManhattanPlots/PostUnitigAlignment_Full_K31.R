# Amy Campbell
# January 2021
# Taking in SAM alignment from BWA Mem of dbgwas-generated unitigs for 221 Isolates (Xanthin Phenotype)
# against SA_502A 

# Required packages
###################
library(tidyverse)
library(dplyr)
library(Rsamtools)
library(qqman)
library(treeio)
library(ggtree)
library(ape)

setwd("~/Desktop/GriceLabGit/DFUStrainsWGS/")

###############
# Read in data 
###############
NodeInfo = data.frame(read_tsv("data/DBGWAS_Xanthin_Full/all_comps_nodes_info.tsv"))
NodeInfoFull = NodeInfo
NodeInfo$Annotations.sep..... = NULL
bamfile = Rsamtools::scanBam("data/DBGWAS_Xanthin_Full/UnitigAlignment.bam")

######################
# QQ plots of p-values 
######################
Unitigs_By_Isolate_Path = "/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/DataInput_DFUStrainWGS_scripts/data/DBGWAS_Xanthin_Full/bugwas_input.all_rows.binary"
Unitigs_By_Isolate_Path_Unique = "/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/DataInput_DFUStrainWGS_scripts/data/DBGWAS_Xanthin_Full/bugwas_input.unique_rows_to_all_rows.binary"

# Checking q-value calculations manually 
# All_Unique_Unitigs = read.table("/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/Results/DBGWAS_Output_Full_31/step1/bugwas_input.unique_rows.binary",header = T )
# All_Unique_Unitigs = All_Unique_Unitigs %>% select(-ps)
# All_Unique_Unitigs = All_Unique_Unitigs %>% rowSums()(
# sum(All_Unique_Unitigs >= 12)
# All_Unique_Unitigs = which(All_Unique_Unitigs %>% filter(>=12)

LineList_UniqueUnitigs = scan(Unitigs_By_Isolate_Path_Unique, what="character", sep="\n")
list_nodes = lapply(LineList_UniqueUnitigs, function(x) paste0("n", (str_split(x," "))[[1]][1]))
rm(LineList_UniqueUnitigs) 
list_nodes = unlist(list_nodes, use.names = FALSE) 
uniques = NodeInfo %>% filter(NodeId %in% list_nodes)



node_info_41  = data.frame(read_tsv("data/all_comps_nodes_info_41.tsv"))
node_info_noPopCorrection = read_tsv("data/all_comps_nodes_info_noPopCorrection.tsv")
# Hacking qqman  qq() fxn to plot qq plots on top of one another for pop-corrected vs. pop-uncorrected
#################################################################################################
ActualFullSetUnitigs = read_delim("/Users/amycampbell/Box/GRICE LAB SHARE/Current lab members/Amy/StaphStrainAnalysis/Results/DBGWAS_Output_Full_31/step2/patterns.txt", " ", col_names=T)
ActualFullSetUnitigs_uncorrected = read_delim("data/patterns_NO_POP_CORRECTION.txt", " ", col_names=T)

pvect_Unfiltered_Corrected = NodeInfoFull$p.value
pvect_Unfiltered_Uncorrected = node_info_noPopCorrection$`p-value`
pvect_BUGWAS = ActualFullSetUnitigs$`p-value`
pvect_BUGWAS_NoCorrection = ActualFullSetUnitigs_uncorrected$`p-value`

pvect_Unfiltered_Corrected <- pvect_Unfiltered_Corrected[!is.na(pvect_Unfiltered_Corrected) & !is.nan(pvect_Unfiltered_Corrected) & !is.null(pvect_Unfiltered_Corrected) & is.finite(pvect_Unfiltered_Corrected) & pvect_Unfiltered_Corrected<1 & pvect_Unfiltered_Corrected>0]
pvect_Unfiltered_Uncorrected <- pvect_Unfiltered_Uncorrected[!is.na(pvect_Unfiltered_Uncorrected) & !is.nan(pvect_Unfiltered_Uncorrected) & !is.null(pvect_Unfiltered_Uncorrected) & is.finite(pvect_Unfiltered_Uncorrected) & pvect_Unfiltered_Uncorrected<1 & pvect_Unfiltered_Uncorrected>0]
pvect_BUGWAS   <- pvect_BUGWAS[!is.na(pvect_BUGWAS) & !is.nan(pvect_BUGWAS) & !is.null(pvect_BUGWAS) & is.finite(pvect_BUGWAS) & pvect_BUGWAS<1 & pvect_BUGWAS>0]
pvect_BUGWAS_NoCorrection <-  pvect_BUGWAS_NoCorrection[!is.na(pvect_BUGWAS_NoCorrection) & !is.nan(pvect_BUGWAS_NoCorrection) & !is.null(pvect_BUGWAS_NoCorrection) & is.finite(pvect_BUGWAS_NoCorrection) & pvect_BUGWAS_NoCorrection<1 & pvect_BUGWAS_NoCorrection>0]


qqobserved_corrected = -log10(sort(pvect_Unfiltered_Corrected,decreasing=FALSE))
qqexpect_corrected = -log10( ppoints(length(pvect_Unfiltered_Corrected) ))
  
  
qqobserved_uncorrected = -log10(sort(pvect_Unfiltered_Uncorrected,decreasing=FALSE))
qqexpect_uncorrected = -log10( ppoints(length(pvect_Unfiltered_Uncorrected) ))

qqobserved_BUGWAS = -log10(sort(pvect_BUGWAS,decreasing=FALSE))
qqexpect_BUGWAS = -log10( ppoints(length(pvect_BUGWAS) ))

qqobserved_BUGWAS_noCorrection = -log10(sort(pvect_BUGWAS_NoCorrection,decreasing=FALSE))
qqoexpect_BUGWAS_noCorrection =  -log10( ppoints(length(pvect_BUGWAS_NoCorrection) ))

# Based on github user "@romunov"'s method for plotting two scatter plots with diff lengths https://gist.github.com/romunov/47f51ccdfe2362c66e60743849fde6b0 
correctedDF = data.frame(expect=qqexpect_corrected, observe=qqobserved_corrected)
uncorrectedDF = data.frame(expect=qqexpect_uncorrected, observe=qqobserved_uncorrected)
bugwasDF = data.frame(expect=qqexpect_BUGWAS, observe=qqobserved_BUGWAS )
bugwasDF_noCorrection = data.frame(expect=qqoexpect_BUGWAS_noCorrection, observe=qqobserved_BUGWAS_noCorrection )


ggplot() + geom_point(data=correctedDF, aes(x=expect, y=observe),  color="black") 
ggplot(correctedDF, aes(x=expect,y=observe)) + geom_point() + geom_point(data=uncorrectedDF, aes(x=expect, y= observe), color="red") +
  geom_point(data=bugwasDF, aes(x=expect, y=observe), color="blue") + 
  geom_abline(slope=1) + annotate("text", label="No correction", color="red", x=1, y=4.8) + annotate("text", label="Correction", color="black", x= 2, y=4) +
  xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") + ggtitle("QQ Plot for GWAS p-values") + annotate("text", label="Actual full set of p-values", color="blue", x=3, y=2)+
  geom_point(data=bugwasDF_noCorrection, aes(x=expect, y=observe), color="green")



qq(ActualFullSetUnitigs$`p-value`)
################################
# Results of association testing
################################
NodeIDs = bamfile[[1]]$qname
positions = bamfile[[1]]$pos
AlignmentFrame = data.frame(rep(0, length(NodeIDs)))
AlignmentFrame$NodeId = NodeIDs   
AlignmentFrame$UnitigPosition= positions

UnalignedUnitigs = AlignmentFrame %>% filter(is.na(UnitigPosition))
UnalignedUnitigsInfo = NodeInfo %>% filter(NodeId %in% UnalignedUnitigs$NodeId)
View(UnalignedUnitigsInfo)

# SA502A genome 
###############
# Chromosome
# 1:2764699
# Plasmid 
# 2764700:2787689

max(AlignmentFrame$UnitigPosition, na.rm=T)
AlignmentFrame$UnitigPosition <- AlignmentFrame$UnitigPosition %>% replace_na(3000000)
AlignmentFrame = AlignmentFrame %>% select(NodeId, UnitigPosition)
GWAS_result = NodeInfo %>% select(NodeId, AlleleFreq, q.Value, p.value, SequenceLength) %>% left_join(AlignmentFrame, by="NodeId")

GWAS_result = GWAS_result %>% filter(!is.na(q.Value))
GWAS_result$logtransformedQ = log10(GWAS_result$q.Value)

unmapped = GWAS_result %>% filter(UnitigPosition == 3000000)
GWAS_result = GWAS_result %>% mutate(ReferenceMap = if_else(UnitigPosition == 3000000, "Unmapped", "Mapped"))

ggplot(GWAS_result, aes(x=UnitigPosition, y=-logtransformedQ)) +
  geom_point(aes(colour=factor(ReferenceMap))) + xlab("Unitig Position on S. aureus 502A") +
  ylab("-Log-10-transformed Q Value") +scale_color_manual(values=c("black", "grey58"))+ geom_line(y=-log10(.01), linetype="dashed", colour="red") +
  geom_line(y=-log10(.05), linetype="dashed", colour="blue") +theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(colour=guide_legend(title="Mapping to SA502A")) + 
  annotate("text", label="q-value = .05", color="blue", x= 80000, y=1.4, hjust=0) +
  annotate("text", label="q-value = .01", color="red", x=80000, y=2.1, hjust=0) +
  xlim(0, 3000000) #+



significant05= (GWAS_result %>% filter(q.Value < .05))

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
  Overlap = (intersect(  ((AnnotationFrame[closest_start_Unitig_start,])$start : (AnnotationFrame[closest_start_Unitig_start,])$end), unitigposition:(unitigposition + unitiglength)))
  if(returnType == 1){
    return((AnnotationFrame[closest_start_Unitig_start,])$attr)
  }else if(returnType == 2){
    return(length(Overlap))
  }else if(returnType==3){
    return((AnnotationFrame[closest_start_Unitig_start,])$start)
  }
}

FindAnnotation(825845,36, SA502A_annotations, 1)  
FindAnnotation(825845,36, SA502A_annotations, 2)  

significant05$Annotation502A = NA
significant05$AnnotationOverlap = NA
significant05$AnnotationLocus = NA

for(r in 1:nrow(significant05)){
  row = significant05[r,]
  significant05[r, "Annotation502A"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 1)
  significant05[r, "AnnotationOverlap"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 2)
  significant05[r, "AnnotationLocus"] <- FindAnnotation(row["UnitigPosition"],row["SequenceLength"], SA502A_annotations, 3)
}
write.csv(significant05, file="QValuesUnder05.csv")
#########################################
# Plot two significant results on a tree
#########################################

# Read in phylogeny
###################
Isolates = "data/DFU_Staph_aureus_isolates.csv"
CoreTree_PosteriorReferences = "data/core_gene_alignment.newick"
ViewTree=treeio::read.jplace("data/core_gene_alignment.jplace")

IsolatesInfo = read.csv(Isolates)
IsolatesInfo$DORN = paste0("DORN", IsolatesInfo$Doern.lab.bank.)
IsolatesInfo$subject_timepoint = paste(IsolatesInfo$patient_id, IsolatesInfo$visit,sep="_")
IsolatesInfo$Genome2=paste0("DORN", IsolatesInfo$Doern.lab.bank.)
IsolatesInfo$Genome1=paste0("DORN", IsolatesInfo$Doern.lab.bank.)
ReferenceNames=c("SA_502A","S_epidermidis",
                 "USA100_AR465", "USA300_FPR3757",
                 "USA400_051", "CC8_NCTC8325",
                 "CC8_Newman", "CC1_MSSA476",
                 "CC1_MW2", "CC22_HO",
                 "CC30_MRSA252", "CC398",
                 "CC398_ATCC6538","CC5_Mu50",
                 "CC5_N315","CC72_CN1", "SA_CFSAN007883", "SA_MCRF184", "SA_AR464", "SA_UP_1150")
PPlacerTree = ggtree::read.tree(CoreTree_PosteriorReferences)
apeRooted = ape::root(PPlacerTree,"S_epidermidis", resolve.root=T)

# Adjust Outgroup length
apeplot <- ggtree::ggtree(apeRooted, size=.25, layout = "rectangular") + geom_tiplab(size=2)
apeplotbackup = apeplot
apeplot$data$DORN = apeplot$data$label
apeplot$data[apeplot$data$label %in% c("S_epidermidis"), "x"] = mean(apeplot$data$x)
apeplot$data$DORN = if_else(apeplot$data$DORN == "CC398_ATCC6538", "ATCC6538", apeplot$data$DORN )
apeplot$data$label = apeplot$data$DORN

UnitigDist = read.table(Unitigs_By_Isolate_Path,header = T )
UnitigDist$node = paste0("n", UnitigDist$ps)


Node_YagU = "n69294"
NodeYagU_distribution = UnitigDist %>% filter(node =="n69294")
TransposedYagU_distribution = data.frame(t(NodeYagU_distribution[,2:(ncol(NodeYagU_distribution) - 1)]))
colnames(TransposedYagU_distribution) = c("n69294")
TransposedYagU_distribution$DORN = row.names(TransposedYagU_distribution)
apeplot$data = apeplot$data %>% left_join(TransposedYagU_distribution, by="DORN")

YagU_Tree = ggtree(apeplot$data) + geom_point(size=5, aes(color=(n69294>0)))+ scale_color_manual(values=c("NA","darkgoldenrod2")) + geom_tiplab() +guides(colour=guide_legend(title="YagU Minor Frequency Variant")) + theme(legend.position=c(.2,.8), legend.title = element_text(size=30), legend.text = element_text(size=25)) #theme(legend.position=c(.2, .8), legend.key.size=unit(20, "cm"))
ggsave(YagU_Tree, file="YagU_Tree.pdf", width=30, height=40)



Node_essG = "n231982"
Node_essG_distribution = UnitigDist %>% filter(node =="n231982")
Transposed_essG_distribution = data.frame(t(Node_essG_distribution[,2:(ncol(Node_essG_distribution) - 1)]))
colnames(Transposed_essG_distribution) = c("n231982")
Transposed_essG_distribution$DORN = row.names(Transposed_essG_distribution)
apeplot$data = apeplot$data %>% left_join(Transposed_essG_distribution, by="DORN")

essG_Tree = ggtree(apeplot$data) + geom_point(size=5, aes(color=(n231982>0)))+ scale_color_manual(values=c("NA","darkgoldenrod2")) + geom_tiplab() +guides(colour=guide_legend(title="EssG Minor Frequency Variant")) + theme(legend.position=c(.2,.8), legend.title = element_text(size=30), legend.text = element_text(size=25)) #theme(legend.position=c(.2, .8), legend.key.size=unit(20, "cm"))
ggsave(essG_Tree, file="EssG_Tree.pdf", width=30, height=40)

ARM_Phenotypes = read.csv("data/Phenotypes_01.15.20.csv")
ARM_Phenotypes$DORN = paste0("DORN", ARM_Phenotypes$DORN)
genotype_yagu_present  = (YagU_Tree$data %>% filter(n69294>0))$DORN

genotype_essG_present  = (essG_Tree$data %>% filter(n231982>0))$DORN
 

ARM_Phenotypes %>% filter(DORN %in% genotype_yagu_present)
ARM_Phenotypes %>% filter(DORN %in% genotype_essG_present)

genotype_essG_present  = (YagU_Tree$data %>% filter(n69294>0))$DORN

Unitigs_By_Isolate_Path
