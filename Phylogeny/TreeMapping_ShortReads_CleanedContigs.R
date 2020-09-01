# Amy Campbell
# 07/2020
# Script to make a tree (ML, core gene alignment-based with the 
# reference genomes/outgroup added a posteriori via pplacer)
# I'm making the full tree labeled by healing outcome, as well as 
# one which is collapsed s.t. I can manually list the very closely related genomes 
# This script also contains the code to output a bar graph of collapsed nodes by count
# of healing outcomes(for each unique subject represented in a given clade), which I'll use
# (in AI) to make a more legible tree of healing outcomes 
# This version is specifically updated to deal with the tree output by Pplacer/RaxML, 
# built on core genome alignments
setwd('/Users/amycampbell/Desktop/Club_Grice/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/')
library("treeio")
library("ggtree")
library("dplyr")
library("ggplot2")

############
# FUNCTIONS
############
SubjTime = function(row){
  if(is.na(row["DORN"])){
    return(NA)
  }else if(is.na(row["patient_id"])){
    return(row["DORN"])
  }else{
    return(paste(row["patient_id"], stringr::str_replace(row["visit"], " ", ""), sep="_"))
  }}

# list of Isolates and what subjects and timepoints they correspond to 
Isolates = "DFU_Staph_aureus_isolates.csv"

# Amelia's phenotypes, which also have the healing outcome/subject for each isolate. 
# Only thing to note here is that 'weeks to heal' is listed as '50' 
# if the wound never healed (which would be NA in the December Gardner metadata)
ARM_Phenotypes = read.csv("Phenotypes_01.15.20.csv")
StaphIsolateDORNs=read.csv("DFU_Staph_aureus_isolates.csv")
CoreTree_PosteriorReferences = "core_gene_alignment.newick"
ViewTree=treeio::read.jplace("core_gene_alignment.jplace")



################################################################
# (1) Assess distances between isolates and reference genomes
#################################################################
# Metadata
###########
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

snpdists_references = read.csv("snp-dists-References.tsv", sep='\t', row.names=1)
rownames(snpdists_references) = colnames(snpdists_references)

snpdists_references$firstgenome = rownames(snpdists_references)
Melted_refs = reshape2::melt(snpdists_references, id.vars=c("firstgenome"))

colnames(Melted_refs) = c("Genome1", "Genome2", "Distance")

# Subset to pairs where the second genome is a reference and the first is not
Melted_refs_Closest = Melted_refs %>% filter(Genome2 %in% ReferenceNames) %>% filter(!(Genome1 %in% ReferenceNames))
# Filter to rows that include the closest (by snp distance) match for each. 
Melted_refs_Closest = Melted_refs_Closest %>% group_by(Genome1) %>% mutate(mindist = min(Distance))
Melted_refs_Closest= Melted_refs_Closest %>% filter(mindist==Distance)

# For those which are equidistant from two closest references, arbitrarily choose the first one in the dataframe. 
Melted_refs_Closest_Summarize = Melted_refs_Closest %>% group_by(Genome1) %>% summarise(Closest=Genome2[1])
Melted_refs_Closest_Summarize = data.frame(Melted_refs_Closest_Summarize)
write.table(Melted_refs_Closest_Summarize,  row.names=F, col.names=F, quote=F, sep="\t", file="ClosestReferences.txt")


####################################################################################
# (2) Plot ML, core genome-SNP-based tree of DORN isolates with S. aureus references
####################################################################################

# IDing leaves on the tree that are within 250 SNPs of each other
# Strictly for purposes of collapsing leaves of the tree 

# Can't be comparing a leaf to itself
Melted_refs_diffgenomes = Melted_refs %>% filter(Genome1 != Genome2)

Melted_refs_diffgenomes = Melted_refs_diffgenomes %>% mutate(AlmostIdentical250= (Distance<250))
AlmostIdenticalPairs_Refs = Melted_refs_diffgenomes %>% filter(AlmostIdentical250)
AlmostIdentical250Pairs_Refs = AlmostIdenticalPairs_Refs
ClusterSizeCounts = AlmostIdentical250Pairs_Refs %>% group_by(Genome1) %>% count(AlmostIdentical250)
AlmostIdentical250RefList = (ClusterSizeCounts %>% arrange(-n))$Genome1

j=1
# Build a set of non-redundant clusters (go through 'genome 1' column arbitrarily, 
# and for each genome that hasn't already been placed in a cluster, make it the representative for a new cluster
# and collapse everything that is <118 snps from it into its cluster

Melted_refs_diffgenomes
AlmostIdentical250ClustersRef = list()
while(length(AlmostIdentical250RefList) > 0){
  RepGenome = AlmostIdentical250RefList[1]
  LittleList=sapply((AlmostIdentical250Pairs_Refs %>% filter(Genome1==RepGenome))$Genome2, function(x) toString(x))
  AlmostIdentical250ClustersRef[[j]] <- append(LittleList, RepGenome)
  j=j+1
  AlmostIdentical250RefList = setdiff(AlmostIdentical250RefList, LittleList)
  AlmostIdentical250RefList = setdiff(AlmostIdentical250RefList, c(RepGenome))
}

length(AlmostIdentical250ClustersRef)

# Read in the tree 
###################
#CoreTree_PosteriorReferences is the tree with the outgroups and ingroups added a posteriori 
PPlacerTree = ggtree::read.tree(CoreTree_PosteriorReferences)
apeRooted = ape::root(PPlacerTree,"S_epidermidis", resolve.root=T)

# Adjust Outgroup length
apeplot <- ggtree(apeRooted, size=.25, layout = "rectangular") + geom_tiplab(size=2)

#S epi branch originally 0.324, becomes 0.0199
# 0.0199/0.324
# Or, 6.15% its original length (approximately 1/16)
# Note this on any publication of the plot 

apeplot$data[apeplot$data$label %in% c("S_epidermidis"), "x"] = mean(apeplot$data$x)
ggsave(apeplot, file="treewithOutScaled_DORNlabels.pdf", width=20, height=20)

# Change labels to be subject, timepoint and save a version where these are the labels 
# as treewithOutScaled_SubjTime.pdf
Plot_Data = data.frame(apeplot$data)

Plot_Data$DORN = Plot_Data$label 
Plot_Data = Plot_Data %>% dplyr::left_join(IsolatesInfo[c("DORN", "patient_id", "visit")], by="DORN")

Plot_Data$TipLABEL = apply(Plot_Data, 1, function(x) SubjTime(x))
Plot_Data$SubjTime= Plot_Data$TipLABEL

# Add healing info
ARM_Phenotypes$patient_id = ARM_Phenotypes$Patient
SubsetVars_Phenotypes =ARM_Phenotypes[c("HealedYorN", "patient_id", "WeekHealed")] %>% distinct()
Plot_Data = Plot_Data %>% left_join(SubsetVars_Phenotypes, by="patient_id")

# Time to heal categories to match LK's
Plot_Data = Plot_Data %>% mutate(HealingTime=case_when(
  (WeekHealed <= 4) ~"<=4 Weeks", 
  ((WeekHealed > 4)  & (WeekHealed <= 8)) ~"<=8 Weeks",
  ((WeekHealed > 8)  & (WeekHealed <= 12)) ~"<=12 Weeks",
  (WeekHealed > 12)  ~">12 Weeks")
)

# Set apeplot data to be this now modified dataframe (all this metadata added)
apeplot$data = Plot_Data

ggsave(apeplot + geom_tippoint(aes(colour=HealingTime)) + scale_colour_manual(values=c("#FFFF99","darkgoldenrod1", "darkorange2", "darkred", "black")), file="AllDORNs_HealingTimeLabels.pdf")

apeplot_subjlabeled = apeplot
apeplot_subjlabeled$data$label = apeplot_subjlabeled$data$SubjTime
ggsave(apeplot_subjlabeled + geom_tippoint(aes(colour=HealingTime)) +  scale_colour_manual(values=c("#FFFF99","darkgoldenrod1", "darkorange2", "darkred", "black")), file="treewithOutScaled_SubjTime.pdf",width=20, height=25)
ggsave(apeplot + geom_tippoint(aes(colour=HealingTime)) +  scale_colour_manual(values=c("#FFFF99","darkgoldenrod1", "darkorange2", "darkred", "black")), file="treewithOutScaled_DORNS.pdf",width=20, height=25)

Selected_ONP = read.table("AllDORNs_ONP.txt", col.names=c("DORN"))
Selected_ONP$Index = row.names(Selected_ONP)

apeplotsave= apeplot
apeplot_data_save = apeplot$data
apeplot_data_save = apeplot_data_save %>% left_join(Selected_ONP, by="DORN")

apeplot_data_save = apeplot_data_save %>% mutate(ONP= !is.na(Index))
apeplotsave$data = apeplot_data_save
ggsave(apeplotsave + geom_tippoint(aes(colour=ONP)) +  scale_colour_manual(values=c("white","green")), file="treewithOutScaled_DORNS_labeledbyONP.pdf",width=20, height=25)
ggsave(apeplot + geom_tippoint(aes(colour=HealingTime)) +  scale_colour_manual(values=c("#FFFF99","darkgoldenrod1", "darkorange2", "darkred", "black")), file="treewithOutScaled_DORNS.pdf",width=20, height=25)

# Identify the two most distant members of each cluster, find their MRCA to be the node at which we collapse that clade
nodelist=c()
for(i in 1:length(AlmostIdentical250ClustersRef)){
  ClusterDF = AlmostIdentical250Pairs_Refs %>% filter(Genome1 %in% AlmostIdentical250ClustersRef[[i]])
  otput = ClusterDF[which.max(ClusterDF$Distance), ]
  collapsenode = phytools::findMRCA(apeRooted, c(toString(otput$Genome1), toString(otput$Genome2)))
  nodelist = c(nodelist, collapsenode)
}

nodelist = unique(nodelist)
nodelist_backup = nodelist
nodelist_compare = nodelist

apeplot_labeled = apeplot

# for each node listed for removal
removelist=c()
for(n in nodelist_backup){
  # get its descendants 
  n_descendants = phytools::getDescendants(as.phylo(apeplot_labeled), node=n)
  for(o in nodelist_compare){
    o_descendants =  phytools::getDescendants(as.phylo(apeplot_labeled), node=o)
    if(o %in% n_descendants){
      removelist = append(removelist, o)
    }
  }}
newnodelist = nodelist[which(!(nodelist %in% removelist))]
    

# Collapse the tree at the nodes listed in "newnodelist"
#########################################################
currentplot=apeplot_labeled
# keep track of children for each collapsed node
ChildrenList = list()
for (i in 1:length(newnodelist)){
  ChildrenList[[i]] = c(phytools::getDescendants(as.phylo(currentplot), node=newnodelist[i]))
  newplot = collapse(currentplot, node=newnodelist[i])
  dataforrearrange = newplot$data
  dataforrearrange = dataforrearrange %>% mutate(isTip = if_else((node==newnodelist[i]), TRUE, isTip)) %>% mutate(label = if_else((node==newnodelist[i]), paste0("NODE", as.character(i)), as.character(label)))
  newplot$data = dataforrearrange
  currentplot=newplot
}
currentplot$data = currentplot$data %>% mutate(Healedby12 =HealingTime %in% c("<=4 Weeks", "<=8 Weeks", "<=12 Weeks"))

# Record children(DORN ID, SubjectTime, and node name) of each node we collapsed the tree at
# Needed for records / to label the collapsed nodes accurately. 
############################################################################################
newdata=(currentplot$data)
newdata$grouping=newdata$DORN
savechildrendata = currentplot$data
for(c in 1:length(ChildrenList)){
  datafr = (savechildrendata %>% filter(node %in% ChildrenList[[c]]))[c("node", "SubjTime","DORN", "HealingTime", "Healedby12")] 
  newdata = newdata %>% mutate(grouping = if_else( (node %in% ChildrenList[[c]]), as.character(c) , grouping))
  write.csv(datafr, file=paste0("NodeChildren_",toString(c)))
}

# Manually ordered in the vertical order by which they appear on currentplot :( 
orderedHeight = c("1", "DORN2178", "2", "DORN1197", "DORN1339", "6", "DORN1334", "11", "12", "4", "DORN2109", "10", "13", "3", "DORN1430", 
                  "DORN1732", "DORN1473", "7", "DORN933", "DORN801", "DORN929", "DORN882", "DORN1679", "DORN900","CC22_HO", "8", "9", "CC398", "DORN1523", "DORN1520", "5", "S_epidermidis")

currentplot$data = currentplot$data %>% mutate(SubjTime = if_else(is.na(SubjTime), label, SubjTime))
currentplot$data$label = currentplot$data$SubjTime
ggsave(ggtree(currentplot$data, size=2) + geom_tiplab(size=15) , file="SimplifiedPlot_labeledBySubjs1.pdf", height=49, width=40)

# Make and save a barplot, ordered the same from top-to-bottom as the tree,
# of healing outcomes for the distribution of subjects represented in each leaf
# of the tree, for overlay onto the tree in adobe illustrator
###########################################################################
newdata$node = NULL
newdata$DORN=NULL

newdata = newdata[c("grouping","patient_id", "Healedby12")]
newdata_melted = newdata %>% filter(grouping!=0) %>% filter(!is.na(patient_id))  %>% reshape2::melt(id.vars=c("grouping","patient_id" ))
BarPlot <- ggplot(newdata_melted, aes(x=as.factor(grouping), fill=value)) + geom_bar() + coord_flip() + scale_fill_manual(values=c("darkred", "#9FE2BF"))  + theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=40), legend.position="bottom", axis.line.x=element_line()) + scale_x_discrete(limits=rev(c("1", "DORN2178", "2", "DORN1197", "DORN1339", "6", "DORN1334", "11", "12", "4", "DORN2109", "10", "13", "3", "DORN1430", 
                                                                                                                                                                                                                                                                                     "DORN1732", "DORN1473", "7", "DORN933", "DORN801", "DORN929", "DORN882", "DORN1679", "DORN900","CC22_HO", "8", "9", "CC398", "DORN1523", "DORN1520", "5", "S_epidermidis")))
BarPlotfinal = BarPlot + theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background  = element_rect(fill = "transparent", colour = NA), axis.ticks.y = element_blank(), axis.text.y=element_blank())
ggsave(BarPlotfinal, file="BarPlotForTreeAnnotation1.pdf", height=15, width=30)



       