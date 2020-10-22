# Amy Campbell
# 07/2020s
# Script to make a tree (ML, core gene alignment-based with the 
# reference genomes/outgroup added a posteriori via pplacer)
# I'm making the full tree labeled by healing outcome, as well as 
# one which is collapsed s.t. I can manually list the very closely related genomes 
# This script also contains the code to output a bar graph of collapsed nodes by count
# of healing outcomes(for each unique subject represented in a given clade), which I'll use
# (in AI) to make a more legible tree of healing outcomes 
# This version is specifically updated to deal with the tree output by Pplacer/RaxML, 
# built on core genome alignments
setwd('/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/')
library("treeio")
library("ggtree")
library("dplyr")
library("ggplot2")
library("dplyr")
library("ape")
library("tibble")

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
Isolates = "data/DFU_Staph_aureus_isolates.csv"

# Amelia's phenotypes, which also have the healing outcome/subject for each isolate. 
# Only thing to note here is that 'weeks to heal' is listed as '50' 
# if the wound never healed (which would be NA in the December Gardner metadata)
ARM_Phenotypes = read.csv("data/Phenotypes_01.15.20.csv")
StaphIsolateDORNs=read.csv("data/DFU_Staph_aureus_isolates.csv")
CoreTree_PosteriorReferences = "data/core_gene_alignment.newick"
ViewTree=treeio::read.jplace("data/core_gene_alignment.jplace")



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

snpdists_references = read.csv("data/snp-dists-References.tsv", sep='\t', row.names=1)
rownames(snpdists_references) = colnames(snpdists_references)

d <- dist(snpdists_references, method = "euclidean")
clustering = hclust(d, method="ward.D")
plot(clustering)
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

# Read in the tree 
###################
#CoreTree_PosteriorReferences is the tree with the outgroups and ingroups added a posteriori 
PPlacerTree = ggtree::read.tree(CoreTree_PosteriorReferences)
apeRooted = ape::root(PPlacerTree,"S_epidermidis", resolve.root=T)

# Adjust Outgroup length
apeplot <- ggtree::ggtree(apeRooted, size=.25, layout = "rectangular") + geom_tiplab(size=2)
apeplotbackup = apeplot

#S epi branch originally 0.324, becomes 0.0199
# 0.0199/0.324
# Or, 6.15% its original length (approximately 1/16)
# Note this on any publication of the plot
apeplot$data$DORN = apeplot$data$label
apeplot$data[apeplot$data$label %in% c("S_epidermidis"), "x"] = mean(apeplot$data$x)
apeplot$data$DORN = if_else(apeplot$data$DORN == "CC398_ATCC6538", "ATCC6538", apeplot$data$DORN )
apeplot$data$label = apeplot$data$DORN
ggsave(apeplot, file="treewithOutScaled_DORNlabels.pdf", width=20, height=20)
# Change labels to be subject, timepoint and save a version where these are the labels
# as treewithOutScaled_SubjTime.pdf
Plot_Data = data.frame(apeplot$data)
Plot_Data = apeplot$data

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

# Collapsing by 2000 snp in a cleaner/more systematic way
##########################################################
apeplot$data = Plot_Data
apeplotbackup = apeplot
apeplotbackup$data$label  = apeplotbackup$data$node
savedata =apeplot$data
phylo_apeplot = as.phylo(apeplotbackup)
ndescendantslist = c()

descendantslist=c()
for(i in 1:nrow(savedata)){
  nodenum = (as.integer(savedata[i, "node"]))
  descendantslist[i] =length(phytools::getDescendants(phylo_apeplot,nodenum))
}

savedata$numDescendants = descendantslist
apeplotbackup$data = savedata

savedata = savedata %>% arrange(-numDescendants)
Melted_refs$Genome1 = if_else(Melted_refs$Genome1=="CC398_ATCC6538", "ATCC6538", as.character(Melted_refs$Genome1))
Melted_refs$Genome2 = if_else(Melted_refs$Genome2=="CC398_ATCC6538", "ATCC6538", as.character(Melted_refs$Genome2))
meltedrefs_diffgenomes = Melted_refs %>% filter(Genome1 != Genome2)
apeplotbackup$data =savedata 

# we'll remove from list_of_parents as we go 
list_of_parents = (savedata %>% arrange(-numDescendants))$node
i=1
collapseapeplot = apeplot 

CollapsedNodeList= c()
ChildrenList2000 = c()
while(length(list_of_parents) > 0){
  nodenum=list_of_parents[1]
  print(nodenum)
  # Base case 
  if(length(list_of_parents) > 1){
    list_of_parents = list_of_parents[2:length(list_of_parents)]
  }else{
    list_of_parents = c()
  }
  # get descendants of the node
  descendants = phytools::getDescendants(as.phylo(collapseapeplot),nodenum)
  descendantGenomes = (collapseapeplot$data %>% filter(node %in% descendants) %>% filter(isTip))$DORN
  #print(descendantGenomes)
  if(nodenum == 302){
    print((meltedsubset$Distance))
  }
  meltedsubset = meltedrefs_diffgenomes %>% filter(Genome1 %in% descendantGenomes) %>% filter(Genome2 %in% descendantGenomes)
  
  # If all leaf descendants of the node are within 1000 SNPs of one another, 

  if(nrow(meltedsubset) > 0){
    if(max(meltedsubset$Distance) < 2000){ 
    print(paste0("NODE ", toString(i)))
    # Save the list of leaves descended from that node
    CollapsedNodeList[[i]] = unique(meltedsubset$Genome1)
    ChildrenList2000[[i]] = descendants
    
    # Collapse the tree at that node
    collapseapeplot = collapse(collapseapeplot, node=nodenum)

    collapseapeplot$data$DORN = if_else(collapseapeplot$data$node == nodenum, paste0("NODE_", toString(i)), collapseapeplot$data$DORN)
    collapseapeplot$data$isTip = if_else(collapseapeplot$data$node == nodenum, TRUE,collapseapeplot$data$isTip)
    # Remove all the descendants from the list of ones to check
    list_of_parents = setdiff(list_of_parents, descendants)
    i = i+1
    }}

  
}
collapseapeplot$data$DORN = if_else(collapseapeplot$data$DORN == "CC398_ATCC6538", "ATCC6538", collapseapeplot$data$DORN )

collapseapeplot$data$label =collapseapeplot$data$DORN
orderedHeight2000 = c("NODE_1", "NODE_2", "NODE_6", "NODE_11","NODE_13", "NODE_4", "DORN2109", "NODE_10", "NODE_12", "NODE_3",
                      "NODE_5", "DORN900", "CC22_HO", "NODE_8", "NODE_9", "CC398", "DORN1523", "DORN1520", "NODE_7", "S_epidermidis"
                      )
collapseapeplot$data$label =collapseapeplot$data$SubjTime
ggsave(ggtree(collapseapeplot$data, size=2) + geom_tiplab(size=15) + xlim(0, .022), file="SimplifiedPlot_labeledBySubjs2000.pdf", height=49, width=40)

# Barplot data 
######################################
barplotdata=(collapseapeplot$data)
barplotdata$grouping1=barplotdata$DORN
barplotdata = barplotdata %>% mutate(HealingTime=case_when(
  (WeekHealed <= 4) ~"<=4 Weeks", 
  ((WeekHealed > 4)  & (WeekHealed <= 8)) ~"<=8 Weeks",
  ((WeekHealed > 8)  & (WeekHealed <= 12)) ~"<=12 Weeks",
  (WeekHealed > 12)  ~">12 Weeks")
)
barplotdata$grouping = barplotdata$grouping1
barplotdata_weekstoheal = barplotdata[c("grouping","patient_id", "HealingTime")]
barplotdatabackup= barplotdata
barplotdata = barplotdata %>% mutate(Healedby12 = HealingTime %in% c("<=4 Weeks", "<=8 Weeks", "<=12 Weeks"))

savechildrendata_barplot = barplotdata

# Iterate the collapsed nodes
for(c in 1:length(ChildrenList2000)){
  print(c)
  datafr = (savechildrendata_barplot %>% filter(node %in% ChildrenList2000[[c]]))[c("node", "SubjTime","DORN", "HealingTime", "Healedby12")] 
  barplotdata$grouping1 = if_else( barplotdata$node %in% ChildrenList2000[[c]], paste0("NODE_", toString(c)), barplotdata$grouping1)
  write.csv(datafr, file=paste0("NodeChildren_Collapsed2000_",toString(c)))
}

barplotdata = barplotdata %>% mutate(SubjTime = if_else(is.na(SubjTime), label, SubjTime))

barplotdata$node = NULL
barplotdata$DORN=NULL
barplotdata$grouping = barplotdata$grouping1

barplotdata
barplotdata[c("grouping","patient_id", "HealedYorN")]
barplotdata_HealedEver = barplotdata[c("grouping","patient_id", "HealedYorN")]
barplotdata_weekstoheal = barplotdata[c("grouping","patient_id", "HealingTime")]
barplotdata = barplotdata[c("grouping","patient_id", "Healedby12")]
barplotdata_melted = barplotdata %>% filter(grouping!=0) %>% filter(!is.na(patient_id))  %>% reshape2::melt(id.vars=c("grouping","patient_id" ))

# Barplot with a count for each isolate
#######################################
BarPlot <- ggplot(barplotdata_melted, aes(x=as.factor(grouping), fill=value)) + geom_bar() + coord_flip() + scale_fill_manual(values=c("darkred", "#9FE2BF"))  + theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=40), legend.position="bottom", axis.line.x=element_line()) + scale_x_discrete(limits=rev(orderedHeight2000))                                              
BarPlotfinal = BarPlot + theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background  = element_rect(fill = "transparent", colour = NA), axis.ticks.y = element_blank(), axis.text.y=element_blank())
ggsave(BarPlotfinal, file="BarPlotForTreeAnnotation_2000_NonUniqueSubj.pdf", height=15, width=30)

BarPlotfinal + theme(axis.text.y=element_text(size=10))
#Take all the clades (after collapsing to <250 SNP distance)


# Barplot with a count for each unique subject represented in a clade
#####################################################################
barplotdata_melted_uniqueSubjects = barplotdata_melted %>% group_by(grouping, patient_id) %>% unique()
BarPlotNew  <- ggplot(barplotdata_melted_uniqueSubjects, aes(x=as.factor(grouping), fill=value)) + geom_bar() + coord_flip() + scale_fill_manual(values=c("darkred", "#9FE2BF"))  + theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=40), legend.position="bottom", axis.line.x=element_line()) + scale_x_discrete(limits=rev(orderedHeight2000))
BarPlotfinalNEW = BarPlotNew + theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background  = element_rect(fill = "transparent", colour = NA), axis.ticks.y = element_blank(), axis.text.y=element_blank())
ggsave(BarPlotfinalNEW, file="BarPlotForTreeAnnotation_UniqueSubjectSet2000.pdf", height=15, width=30)

# With Healing time (4 week increments)
barplotdata_weekstoheal_melt = barplotdata_weekstoheal %>% filter(grouping!=0) %>% filter(!is.na(patient_id))  %>% reshape2::melt(id.vars=c("grouping","patient_id"))
barplotdata_melted_uniqueSubjects_HealingTime = barplotdata_weekstoheal_melt %>% group_by(grouping, patient_id) %>% unique()

BarPlotWeeksToHeal <- ggplot(barplotdata_melted_uniqueSubjects_HealingTime, aes(x=as.factor(grouping), fill=value )) + geom_bar() + coord_flip() + scale_fill_manual(values=(c("#FFFF99","darkgoldenrod1", "darkorange2", "darkred"))) + theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=40), legend.position="bottom", axis.line.x=element_line()) + scale_x_discrete(limits=rev(orderedHeight2000))
BarPlotWeeksToHeal$data$value =  factor(BarPlotWeeksToHeal$data$value, levels=(c("<=4 Weeks", "<=8 Weeks","<=12 Weeks", ">12 Weeks")))

ggsave(BarPlotWeeksToHeal +theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background  = element_rect(fill = "transparent", colour = NA), axis.ticks.y = element_blank(), axis.text.y=element_blank()), file="BarPlotForTreeAnnotation_UniqueSubjectSet2000_WeeksHealed.pdf", height=15, width=30)


# Staphyloxanthin box plots 
###############################
ARM_Phenotypes_Xanthins = ARM_Phenotypes[c("DORN", "XanthinPhenotype")]
orderedHeight2000
ChildrenList2000
CollapsedNodeList

# building staphyloxanthin dataframe
phenotypeDF = data.frame()
for(h in orderedHeight2000){
  print(h)
  if(grepl("DORN",h)){
    PhenotypeRows = ARM_Phenotypes_Xanthins %>% filter(DORN==h)
    PhenotypeRows$grouping = h
    phenotypeDF = rbind(phenotypeDF, PhenotypeRows)
  }else{
    if(grepl("NODE_",h)){
      h_integer = as.integer(stringr::str_replace(h, "NODE_", ""))
      nodechildlist= CollapsedNodeList[[h_integer]]
      PhenotypeRows = ARM_Phenotypes_Xanthins %>% filter(DORN %in% nodechildlist)
      PhenotypeRows$grouping = h 
      phenotypeDF = rbind(phenotypeDF, PhenotypeRows)
    }
  }
}
Xanthin_boxplot  <- ggplot(phenotypeDF, aes(x=as.factor(grouping), y=XanthinPhenotype)) + geom_boxplot( fill="yellow3", color="yellow4", alpha=.8) + coord_flip() + theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=40), legend.position="bottom", axis.line.x=element_line()) + scale_x_discrete(limits=rev(orderedHeight2000))
ggsave(Xanthin_boxplot +theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background  = element_rect(fill = "transparent", colour = NA), axis.ticks.y = element_blank(), axis.text.y=element_blank()), file="BarPlotForTreeAnnotation_UniqueSubjectSet2000_Xanthin.pdf", height=15, width=30)






