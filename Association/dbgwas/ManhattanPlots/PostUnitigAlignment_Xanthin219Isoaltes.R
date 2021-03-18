# Amy Campbell
# January 2021
# Taking in SAM alignment from BWA Mem of dbgwas-generated unitigs against SA_502A 

# Required packages
###################
library(tidyverse)
library(dplyr)
library(Rsamtools)
library(qqman)
library("treeio")
library("ggtree")
library("ape")
# Read in DBGWAS Tabular Output
###############################
NodeInfo = data.frame(read_tsv("data/all_comps_nodes_info.tsv"))
NodeInfoFull = NodeInfo
NodeInfo$Annotations.sep..... = NULL
bamfile = Rsamtools::scanBam("data/UnitigAlignment.bam")

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
GWAS_resultFilteredAlleleFreq = GWAS_result %>% filter(AlleleFreq < 208 & AlleleFreq > 11)
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



# Phylogeny
############
Nodes_EssG_Signif = c("n304326", "n231970")

Unitigs_By_Isolate_Path = "/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/Results/DBGWAS_Output_WithAnnotations/step1/bugwas_input.all_rows.binary"
UnitigDist = read.table(Unitigs_By_Isolate_Path,header = T )
UnitigDist$node = paste0("n", UnitigDist$ps)

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

# essG
###############
# Are the two essG-mapped or annotated unitigs identically distributed? 
unique(colSums(UnitigDist %>% filter( node %in% Nodes_EssG_Signif) %>% select(-node, -ps)))
# Yes. Therefore, we can use one of them (say, n231970) to mark their presence on the tree. 
Node231970_distribution = UnitigDist %>% filter (node =="n231970")
TransposedNode231970_distribution = data.frame(t(Node231970_distribution[,2:(ncol(Node231970_distribution) - 1)]))

colnames(TransposedNode231970_distribution) = c("n231970")
TransposedNode231970_distribution$DORN = row.names(TransposedNode231970_distribution)
apeplot$data = apeplot$data %>% left_join(TransposedNode231970_distribution, by="DORN")
Node231970_distribution_Tree = ggtree(apeplot$data) + geom_point(aes(color=(n231970>0),shape=(n231970>0), size=10))+ scale_color_manual(values=c("NA","darkgoldenrod2")) + geom_tiplab()
ggsave(Node231970_distribution_Tree, width=40,height=40, file="Node231970_distribution_Tree.pdf")

# yagU
Node69289dist = UnitigDist %>% filter(node=="n69289")
TransposedNode69289dist = data.frame(t(Node69289dist[,2:(ncol(Node69289dist) - 1)]))
colnames(TransposedNode69289dist) = c("n69289")
TransposedNode69289dist$DORN = row.names(TransposedNode69289dist)

apeplot$data = apeplot$data %>% left_join(TransposedNode69289dist, by="DORN")
TransposedNode69289dist_Tree = ggtree(apeplot$data) + geom_point(aes(color=(n69289>0),shape=(n69289>0), size=10))+ scale_color_manual(values=c("NA","lightgoldenrod1")) + geom_tiplab()
ggsave(TransposedNode69289dist_Tree, width=40,height=40, file="Node69289_yagU_distribution_Tree.pdf")
TransposedNode69289dist_Tree


Node69289dist = UnitigDist %>% filter(node=="n69289")

Node69289dist
max(UnitigDist %>% select(-node, -ps ) %>% rowSums())



UnitigDist %>% distinct_at(vars(colnames(UnitigDist %>% select(-node, -ps))))

