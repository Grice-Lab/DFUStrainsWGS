# amy campbell 
# Phylogeny of DORN isolates + their references (247 genomes total)
# Updated 05-2022
setwd('/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/')
# library("treeio")
# library("dplyr")
# library("ggplot2")
# library("dplyr")
# library("ape")
# library("tibble")


library("ape")
library("ggtree")
library("dplyr")
library("ggplot2")
library("ggtreeExtra")

TreeFilePath = "data/Phylogeny2022Data/RAxML_bestTree.RaxMLTree2022.newick"
clonalframetreepath="data/Phylogeny2022Data/clonalframe_tree.newick.labelled_tree.newick"
CCmapping="data/Phylogeny2022Data/CCMapPlotting.csv"

StaphIsolateDORNs=read.csv("data/DFU_Staph_aureus_isolates.csv")
UpdatedPhenotypes=read.csv("data/staphyloxanthin_paper_data.csv")
UpdatedPhenotypes$DORN=paste0("DORN", UpdatedPhenotypes$DORN)

StaphIsolateDORNs$DORN=paste0("DORN", StaphIsolateDORNs$Doern.lab.bank.)
CCgroupings=read.csv(CCmapping)
# Read in the clonalframeML-corrected tree :) 
###############################################
TreeObjectCF = ggtree::read.tree(clonalframetreepath)



# + scale_fill_manual(values=c("darkred", "#9FE2BF"))

# Root at CC398
###############
apeRootedCF= ape::root(TreeObjectCF,"CC398", resolve.root=T)

apeplotCF <- ggtree::ggtree(apeRootedCF, size=.25, layout = "rectangular") + geom_tiplab(size=2)
apeplotcirc<- ggtree::ggtree(apeRootedCF, size=.25, layout = "circular") + geom_tiplab(size=4.5)

AllRefs = sort((apeplotcirc$data %>% filter(isTip==TRUE))$label)
AllRefs[!grepl(AllRefs, pattern="DORN")]


circdata = apeplotcirc$data
circdata$DORN = circdata$label
circdata = circdata %>% left_join(CCgroupings, by="DORN")

circdata = circdata %>% mutate(CCRefs =  case_when(label=="CC1_MSSA476" ~"CC1",
                                                           label=="CC1_MW2" ~ "CC1",
                                                           label=="CC1-ST188" ~ "CC1", 
                                                           label=="CC12" ~ "CC12", 
                                                           label=="CC133"~"CC133",
                                                           label=="CC15" ~ "CC15", 
                                                           label=="CC22_HO" ~"CC22", 
                                                           label=="CC30_MRSA252" ~ "CC30", 
                                                           label=="CC398" ~ "CC398", 
                                                           label=="CC45" ~ "CC45", 
                                                           label=="CC5_Mu50"~ "CC5",
                                                           label=="CC5_N315"~ "CC5",
                                                           label=="CC59" ~ "CC59", 
                                                           label=="CC7" ~ "CC7", 
                                                           label=="CC72_CN1" ~ "CC72", 
                                                           label=="CC8_NCTC8325" ~ "CC8",
                                                           label=="CC8_Newman" ~ "CC8",
                                                           label=="CC80" ~ "CC80",
                                                           label=="CC9" ~ "CC9", 
                                                           label== "CC97_ATCC6538" ~ "CC97",
                                                           label== "CC97_Newbould305" ~ "CC97", 
                                                           label=="SA_502A" ~ "CC5",
                                                           label=="SA_AR464" ~ "CC45",
                                                           label=="ST20"  ~ "CC20",
                                                           label=="USA100_AR465"~ "CC5",
                                                           label=="USA300_FPR3757" ~ "CC8",
                                                           label=="USA400_051" ~ "CC1",
                                                            TRUE ~ as.character(CCLabel)))
randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#F6BE00",
                "#33A02C","#E6F5C9")
circdataNullLabel=circdata
circdataNullLabel$label=""
circdataNullLabel$TipLABEL=NULL
CircularPlotColors = ggtree(circdataNullLabel, layout="circular",size=.5) + geom_fruit(aes(y=CCRefs,x=0,fill=CCRefs), geom=geom_tile,offset=5, size=.1)  + geom_tiplab(size=0)
CircularPlotColors = CircularPlotColors + scale_fill_manual(values=randpalette18)

apeplotcircForPlot <- ggtree::ggtree(apeRootedCF, size=.25, layout = "circular") + geom_tiplab(size=4.5)

ggsave(CircularPlotColors, file="ColorsCCplotUTD.pdf", width=26, height=25)
ggsave(apeplotcirc, file="CircularPlotCF_UTD.pdf", width=25, height=25)


# Patient 149 
Stx = UpdatedPhenotypes %>% select(staphyloxanthin,DORN)



apeplotcircXanthin<- ggtree::ggtree(apeRootedCF, size=.25, layout = "circular", branch.length = "none") + geom_tiplab(size=1.7)

apeplotcircXanthin$data$DORN = apeplotcirc$data$label
apeplotcircXanthin$data = apeplotcircXanthin$data %>% left_join(CCgroupings, by="DORN")
apeplotcircXanthin$data  = apeplotcircXanthin$data %>% mutate(CCRefs =  case_when(label=="CC1_MSSA476" ~"CC1",
                                                               label=="CC1_MW2" ~ "CC1",
                                                               label=="CC1-ST188" ~ "CC1", 
                                                               label=="CC12" ~ "CC12", 
                                                               label=="CC133"~"CC133",
                                                               label=="CC15" ~ "CC15", 
                                                               label=="CC22_HO" ~"CC22", 
                                                               label=="CC30_MRSA252" ~ "CC30", 
                                                               label=="CC398" ~ "CC398", 
                                                               label=="CC45" ~ "CC45", 
                                                               label=="CC5_Mu50"~ "CC5",
                                                               label=="CC5_N315"~ "CC5",
                                                               label=="CC59" ~ "CC59", 
                                                               label=="CC7" ~ "CC7", 
                                                               label=="CC72_CN1" ~ "CC72", 
                                                               label=="CC8_NCTC8325" ~ "CC8",
                                                               label=="CC8_Newman" ~ "CC8",
                                                               label=="CC80" ~ "CC80",
                                                               label=="CC9" ~ "CC9", 
                                                               label== "CC97_ATCC6538" ~ "CC97",
                                                               label== "CC97_Newbould305" ~ "CC97", 
                                                               label=="SA_502A" ~ "CC5",
                                                               label=="SA_AR464" ~ "CC45",
                                                               label=="ST20"  ~ "CC20",
                                                               label=="USA100_AR465"~ "CC5",
                                                               label=="USA300_FPR3757" ~ "CC8",
                                                               label=="USA400_051" ~ "CC1",
                                                               TRUE ~ as.character(CCLabel)))




apeplotcircXanthin$data = apeplotcircXanthin$data %>% left_join(Stx, by="DORN")
apeplotcircXanthin$data = apeplotcircXanthin$data %>% left_join(StaphIsolateDORNs, by="DORN")
apeplotcircXanthin$data$label=apeplotcircXanthin$data$patient_id

# Plot with bars for staphyloxanthin production & CC 
#####################################################
WithXanthin = apeplotcircXanthin + geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=4, pwidth=6, color=NA) +
  scale_fill_manual(values=randpalette18)+  new_scale_fill()+
  geom_fruit(aes(x=(staphyloxanthin)),fill="#B8860B", geom=geom_bar, stat="identity", orientation="y",offset=.1, pwidth=20)

ggsave(WithXanthin, file="TestXanthinTree.pdf", width=11,height=11)


circdata_healing = ggtree(circdata, layout="circular", size=.5) + geom_tiplab(size=0) +
  geom_tippoint(aes(color = HealedBy12),  size=2) + scale_color_manual(values=c("darkred", "#9FE2BF"),na.value ="#00000000" ) +
  new_scale_fill()+ geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=8,color=NA, offset=.05) + scale_fill_manual(values=randpalette18)

circdata_healing


ggsave(circdata_healing, file="CircularPlotCF_HealingAlone.pdf", width=10, height=10)
ggsave(circdata_healing, file="CircularPlotCF_Healing.pdf", width=12, height=12)

# Make simplified version of the plot with healing outcomes
############################################################
UpdatedPhenotypes = UpdatedPhenotypes %>% mutate(HealedBy12 = if_else(week_healed > 12 | is.na(week_healed), "No", "Yes"))
HealingOutcomes = UpdatedPhenotypes %>% select(DORN, HealedBy12,patient)
circdata = circdata %>% left_join(HealingOutcomes, by="DORN")


# Choose nodes to collapse 
apeplotBar = ggtree::ggtree(circdata, size=.5, layout = "rectangular") + geom_tiplab() 
apeplotBar$data$label = apeplotBar$data$CCRefs
apeplotBarForCollapse = apeplotBar  + geom_text2(aes(subset=!isTip, label=node)) +  geom_tiplab(size=2.5)
ggsave(apeplotBarForCollapse, file="CClabeledNodesForCollapse.pdf", height=30, width=30)

# Nodes to collapse
# 258 -- CC8
# 298 -- CC15
# 311 -- CC20
# 253 -- CC1
# 462 -- CC97
# 323 -- CC7
# 322 -- CC12
# 446 --  CC72
# 458 -- CC9
# 355 -- CC5
# 348 -- CC59
# 330 -- CC45
# 334 -- CC30
NodesCollapse = c(258,298,311, 253, 462, 323, 322, 446, 458, 355, 348, 330, 334)


# pasted:

phyloObjectPlot = as.phylo(apeplotBar)

apeplotBardata = data.frame(apeplotBar$data)
apeplotBardata$NodeNum = "NA"

descendantslist=c()
for(nodenum in NodesCollapse){
  apeplotBardata = apeplotBardata %>% mutate(NodeNum = if_else(node %in%  append(phytools::getDescendants(phyloObjectPlot,nodenum), nodenum), paste0("NODE_", nodenum), as.character(NodeNum)))
  
}
#fill the rest in
apeplotBardata = apeplotBardata %>% mutate(NodeNum = if_else(NodeNum=="NA", paste0("NODE_", node), NodeNum))
collapsedTreeTEST = collapse(apeplotBar, node=443)

apeplotBar$data$NodeNum= apeplotBardata$NodeNum

collapsedTree = apeplotBar
for(nodenum in NodesCollapse){
  print(nodenum)
  collapsedTree = collapse(collapsedTree, node=nodenum)
  collapsedTree$data$isTip[collapsedTree$data$node==nodenum] <- TRUE
}

collapsedTree$data$label = collapsedTree$data$NodeNum
ggtree(collapsedTree$data)


NodeOrder = (collapsedTree$data %>% filter(isTip==TRUE) %>% filter(!is.na(y)) %>% arrange(-y))$NodeNum


barplotdata = collapsedTree$data
barplotdata$grouping=barplotdata$NodeNum

barplotdata = barplotdata[c("grouping","patient", "HealedBy12")]
barplotdata_melted = barplotdata %>% filter(grouping!=0) %>% filter(!is.na(patient))  %>% reshape2::melt(id.vars=c("grouping","patient" ))

barplotdata_melted_uniqueSubjects = barplotdata_melted %>% group_by(grouping, patient) %>% unique()


BarPlotNew  <- ggplot(barplotdata_melted_uniqueSubjects, aes(x=as.factor(grouping), fill=(value))) + geom_bar() + coord_flip() + scale_fill_manual(values=c("darkred", "#9FE2BF"))  + theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=40), legend.position="bottom", axis.line.x=element_line())+ scale_x_discrete(limits=rev(NodeOrder)) 
BarPlotPrint = BarPlotNew + theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background  = element_rect(fill = "transparent", colour = NA), axis.ticks.y = element_blank(), axis.text.y=element_blank())

BarPlotPrint

collapsedTree = ggtree(collapsedTree$data, size=1)
ggsave(file="CollapsedSBDRCTree_labelednodeUTD5-22.pdf", collapsedTree +geom_treescale(linesize = 4) + geom_tiplab(), width=15, height=20) 

ggsave(BarPlotPrint, file="BarplotforHealingByCC.pdf")


#
