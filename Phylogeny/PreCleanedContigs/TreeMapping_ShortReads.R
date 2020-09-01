# Amy Campbell
# 03/2020
# Script to make a tree (ML, core gene alignment-based with the 
# reference genomes/outgroup added a posteriori via pplacer)
# I'm making the full tree labeled by healing outcome, as well as 
# one which is collapsed s.t. I can manually list the very closely related genomes 
# This script also contains the code to output a bar graph of collapsed nodes by count
# of healing outcomes(for each unique subject represented in a given clade), which I'll use
# (in AI) to make a more legible tree of healing outcomes 

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

# Tree built *with* references as part of the core alignment
CoreTreeML = "/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/RAxML_bestTree.References_align_Tree.newick"

# list of Isolates and what subjects and timepoints they correspond to 
Isolates = "/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/DFU_Staph_aureus_isolates.csv"
CoreTree_PosteriorReferences = "/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/PplacerOutputTOG.newick"

ARM_Phenotypes = read.csv("/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/Phenotypes_01.15.20.csv")

LK_metadata=read.csv("/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/metamapLK.csv")


IsolatesInfo = read.csv(Isolates)
IsolatesInfo$DORN = paste0("DORN", IsolatesInfo$Doern.lab.bank.)
ViewTree=treeio::read.jplace("/Users/amycampbell/Desktop/DFU/core_gene_alignment_UpperCase.jplace")

#CoreTree_PosteriorReferences is the tree with the outgroups and ingroups added a posteriori 
PPlacerTree = ggtree::read.tree(CoreTree_PosteriorReferences)
apeRooted = ape::root(PPlacerTree,"S_EPIDERMIDIS", resolve.root=T)

# S. epidermidis outgroup branch length is huge. Figure out what branch # is and convert it 
# to be the average of all the other lengths
plottree = ggtree(apeRooted) + geom_tiplab()
edge=data.frame()
edgeframe=data.frame(apeRooted$edge, edge_num=1:nrow(apeRooted$edge))
colnames(edgeframe) = c("parent", "node", "edge_num")
plottree %<+% edgeframe + geom_label(aes(x=branch, label=edge_num))

# It's #470
# Adjust Outgroup length
apeplot <- ggtree(apeRooted, size=.25, layout = "rectangular") + geom_tiplab(size=2)

apeplot$data[apeplot$data$label %in% c("S_EPIDERMIDIS"), "x"] = mean(apeplot$data$x)
ggsave(apeplot, file="treewithOutScaled_DORNlabels.pdf", width=20, height=20)

# Change labels to be subject, timepoint
Plot_Data = data.frame(apeplot$data)
Plot_Data$DORN = Plot_Data$label 
Plot_Data = Plot_Data %>% dplyr::left_join(IsolatesInfo[c("DORN", "patient_id", "visit")])
Plot_Data$DORN
Plot_Data$TipLABEL = apply(Plot_Data, 1, function(x) SubjTime(x))
apeplot$data$label = Plot_Data$TipLABEL
apeplot$data$DORN = Plot_Data$DORN
ggsave(apeplot, file="treewithOutScaled_SubjTime.pdf",width=20, height=25)

# Annotate with healing info 
ARM_Phenotypes$patient_id  = ARM_Phenotypes$Patient
SubsetVars_Phenotypes =ARM_Phenotypes[c("HealedYorN", "patient_id", "WeekHealed")]
Plot_Data_Updated = Plot_Data %>% left_join(SubsetVars_Phenotypes, by="patient_id")


apeplot_labeledheal <- apeplot

Plot_Data_Updated = Plot_Data_Updated %>% mutate(HealingTime=case_when(
     (WeekHealed <= 4) ~"<=4 Weeks", 
     ((WeekHealed > 4)  & (WeekHealed <= 8)) ~"<=8 Weeks",
     ((WeekHealed > 8)  & (WeekHealed <= 12)) ~"<=12 Weeks",
     (WeekHealed > 12)  ~">12 Weeks")
   )
apeplot_labeledheal$data = Plot_Data_Updated
levels(as.factor(apeplot_labeledheal$data$HealingTime))
apeplot_labeledheal$data$label = apeplot_labeledheal$data$TipLABEL
apeplot_labeledheal$data$HealingTime <- factor(apeplot_labeledheal$data$HealingTime, levels=c("<=4 Weeks", "<=8 Weeks", "<=12 Weeks", ">12 Weeks"))

apeplot_labeledheal = apeplot_labeledheal + geom_tippoint(aes(colour=HealingTime)) + scale_colour_manual(values=c("#FFFF99","darkgoldenrod1", "darkorange2", "darkred", "black"))
ggsave(apeplot_labeledheal, file="treewithOutScaled_LabelHeal.pdf", width=20, height=20)

# Collapse busy clades
######################
# Nodes to 'collapse' and manually label instead
ggtree(apeRooted) + geom_label(aes(label=node))
apeRootedCollapsedALLbut = ape::drop.tip(apeRootedCollapsedALLbut, "S_EPIDERMIDIS")

apeRootedCollapsedALLbut=(apeplot %>% collapse(node=289)) #+geom_label(aes(label=node))
View((apeplot$data)[c("node", "label")])


ggsave(apeRootedCollapsedALLbut, height=30,width=40, file="test_label.pdf")
sort(phytools::getDescendants(apeRooted, node=294))
apeRooted = drop

phytools::findMRCA(apeRooted, c("DORN1410", "DORN568"))
phytools::findMRCA(apeRooted, c("DORN1520", "DORN568"))
phytools::findMRCA(apeRooted, c("DORN1520", "S_EPIDERMIDIS"))
phytools::findMRCA(apeRooted, c("DORN1520", "CC30_MRSA252"))

apeRooted$tip.label

clade327 = groupClade(apeplot, .node=327)
View(clade327)
apeplot <- collapse(apeplot, node=294)
clade327$data$label


apeplot$data[apeplot$data$label %in% c("S_EPIDERMIDIS"), "x"] = mean(apeplot$data$x)
ggsave(apeplot, file="treewithOutScaled_DORNlabels.pdf", width=20, height=20)


apeplot0 <- ggtree(apeRooted, size=.5, layout = "rectangular") + geom_tiplab(size=11.5)
apeplot0$data[apeplot$data$label %in% c("S_EPIDERMIDIS"), "x"] = mean(apeplot0$data$x)

Plot_Data = data.frame(apeplot0$data)
Plot_Data$DORN = Plot_Data$label
Plot_Data = Plot_Data %>% dplyr::left_join(IsolatesInfo[c("DORN", "patient_id", "visit")])
Plot_Data$DORN
Plot_Data$TipLABEL = apply(Plot_Data, 1, function(x) SubjTime(x))
apeplot0$data$label = Plot_Data$TipLABEL
apeplot0$data$DORN = Plot_Data$DORN


apeplot0
apeplot1 <- collapse(apeplot0, node=327)
apeplot2 <- collapse(apeplot1, node=424)
apeplot3 <- collapse(apeplot2, node=413)
apeplot4 <- collapse(apeplot3, node=322)
apeplot5 <- collapse(apeplot4, node=302)
apeplot6 <- collapse(apeplot5, node=261)
apeplot7 <- collapse(apeplot6, node=269)
apeplot8 <- collapse(apeplot7, node=279)
apeplot9 <- collapse(apeplot8, node=289)
apeplot10 <- collapse(apeplot9, node=315)
apeplot11 <- collapse(apeplot10, node=284)
apeplot12 <- collapse(apeplot11, node=318)
ggsave(apeplot12 , height=40, width=49.5, file="CollapsedTree.pdf")

CollapsedNodes = c(327, 424, 413, 322, 302, 261, 269, 279, 289, 315, 284, 318)
bigdf = data.frame()
for(c in CollapsedNodes){
  print(c)
  descendants = phytools::getDescendants(apeRooted, node=c)
  rows = apeplot0$data %>% filter(node %in% descendants)
  all_timeSubj = (unique(unlist(rows["label"])))
  littledf = data.frame((1:length(all_timeSubj))) 
  littledf$timesubj = (all_timeSubj)
  littledf$time  = sapply(littledf$timesubj, function(x) (strsplit(x, "_")[[1]][2]))
  littledf$subject  = sapply(littledf$timesubj, function(x) (strsplit(x, "_")[[1]][1]))
  littledf$healed = sapply(littledf$subject, function(x) toString(unique((ARM_Phenotypes[which(ARM_Phenotypes$Patient==x), "HealedYorN"]))))
  littledf$node = rep(c, length(all_timeSubj))
  bigdf = rbind(bigdf, littledf)

  write.csv(littledf, paste0( toString(c), "_Children.csv"))
  
}

biggerdf = bigdf %>% group_by(node, subject) %>% summarize(max= max(time), healed=first(healed))
countdf = biggerdf %>% group_by(node, healed) %>% summarize(healing = n() )
biggerdf = biggerdf %>% filter(healed!="")
melted = reshape2::melt(biggerdf, id.vars=c("node", "max", "subject"))


BarPlot <- ggplot(melted, aes(x=as.factor(node), fill=value)) + geom_bar() + coord_flip() + scale_fill_manual(values=c("darkred", "#FFFF99"))

colnames(countdf) = c("Node", "Healed", "Count")
apeplot_dorns <- apeplot
apeplot_dorns$data$label = Plot_Data$DORN
p <- plot(apeplot)
ape::edgelabels()
collapse(apeplot, node=326)




ggsave((apeplot + geom_label(aes(label=node))), width=45, height=20,file="labelednodes.pdf")




       