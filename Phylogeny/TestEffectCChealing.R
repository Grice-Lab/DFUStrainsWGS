# Amy Campbell
# 221 isolates tree from alignment including S. epidermidis and all 
# for SBDRC retreat 2022
# (PGAP annotated) 

setwd('/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/')

library("ggtree")
library("dplyr")
library("ggplot2")
library("dplyr")
library("ape")
library("tibble")
library("ggtreeExtra")


# Remove DORN1197, DORN1471, DORN575, DORN22

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


CCestimates = read.csv("mappings/PreliminaryClassificationsCCsSBDRC.csv") 
CCestimates = CCestimates %>% select(Genome, CC)
CCestimates$DORN = CCestimates$Genome
Isolates = "data/DFU_Staph_aureus_isolates.csv"
ARM_Phenotypes = read.csv("data/Phenotypes_01.15.20.csv")
StaphIsolateDORNs=read.csv("data/DFU_Staph_aureus_isolates.csv")
UpdatedPhenotypes=read.csv("data/staphyloxanthin_paper_data.csv")





# Note: These were based on 20 random starts as opposed to 100: 
CoreTree_PosteriorReferences = "data/207Isolates/PPlacerTree207Isolates20Reps.newick"
SBDRCTree = "data/RAxML_bestTree.RaxMLTreeWithRefsSBDRCpreliminary.newick"

# this one is based on the 100 random starts
# CoreTree_PosteriorReferences = "data/207Isolates/core_gene_alignment207.newick"

#######################################
# Set up tree and map metadata onto it
#######################################

# Some metadata
IsolatesInfo = read.csv(Isolates)
IsolatesInfo$DORN = paste0("DORN", IsolatesInfo$Doern.lab.bank.)
IsolatesInfo$subject_timepoint = paste(IsolatesInfo$patient_id, IsolatesInfo$visit,sep="_")
IsolatesInfo$Genome2=paste0("DORN", IsolatesInfo$Doern.lab.bank.)
IsolatesInfo$Genome1=paste0("DORN", IsolatesInfo$Doern.lab.bank.)


IsolatesInfo = IsolatesInfo %>% left_join(CCestimates, by="DORN")

IsolatesInfoSubset = IsolatesInfo %>% select(patient_id,DORN,visit, CC) %>% tidyr::drop_na()


# # note CC398_ATCC6538 is actually not CC398 -- it needs to always be relabeled to ATCC6538 (it's Rosenbach)
# ReferenceNames=c("USA100_AR465", "CC72_CN1", "CC1_MW2", "USA300_FPR3757","CC15",
#                  "CC59", "CC398", "SA_CFSAN007883s", "CC30_MRSA252", "CC5_Mu50",
#                  "CC8_NCTC8325", "CC1_MSSA476", "SA_MCRF184", "SA_AR464", "CC398_ATCC6538",
#                  "CC12", "USA400_051", "CC8_Newman", "SA_502A", "CC22_HO", "S_epidermidis",
#                  "CC80", "CC5_N315", "CC133", "SA_UP_1150")

ReferenceNames = c("CC12", "CC133","CC15", "CC1_MSSA476","CC1_MW2", "CC1-ST188", "CC22_HO",
  "CC30_MRSA252", "CC398", "CC45", "CC59", "CC5_Mu50", "CC5_N315","CC72_CN1",
  "CC7", "CC80", "CC8_NCTC8325", "CC8_Newman", "CC97_ATCC6538", "CC97_Newbould305",
  "SA_502A", "SA_AR464", "S_epidermidis", "ST20", "USA100_AR465", "USA300_FPR3757", "USA400_051")


PPlacerTree = ggtree::read.tree(SBDRCTree)
apeRooted = ape::root(PPlacerTree,"S_epidermidis", resolve.root=T)


apeplot <- ggtree::ggtree(apeRooted, size=.15, layout = "rectangular") + geom_tiplab(size=2)
apeplotbackup = apeplot
apeplot$data$DORN = apeplot$data$label
apeplot$data[apeplot$data$label %in% c("S_epidermidis"), "x"] = mean(apeplot$data$x)

ggsave(apeplot, file="VeryPreliminarySBDRC.pdf", width=25, height=25)
# Change mislabeled ATCC6538
apeplot$data$DORN = if_else(apeplot$data$DORN == "CC398_ATCC6538", "ATCC6538", apeplot$data$DORN )
apeplot$data$label = apeplot$data$DORN


Plot_Data = data.frame(apeplot$data)
Plot_Data = apeplot$data

Plot_Data$DORN = Plot_Data$label
Plot_Data = Plot_Data %>% dplyr::left_join(IsolatesInfoSubset[c("DORN", "patient_id", "visit", "CC")], by="DORN")

Plot_Data$TipLABEL = apply(Plot_Data, 1, function(x) SubjTime(x))
Plot_Data$SubjTime= Plot_Data$TipLABEL

# to prune prior to plotting for SBDRC retreat
RemoveTips  = setdiff(Plot_Data$DORN, append(IsolatesInfoSubset$DORN, "S_epidermidis"))

# Add healing info
ARM_Phenotypes$patient_id = ARM_Phenotypes$Patient
SubsetVars_Phenotypes =ARM_Phenotypes[c("HealedYorN", "patient_id", "WeekHealed")] %>% distinct()
Plot_Data = Plot_Data %>% left_join(SubsetVars_Phenotypes, by="patient_id")
Plot_Data = Plot_Data %>% mutate(HealingTime=case_when(
  (WeekHealed <= 4) ~"<=4 Weeks",
  ((WeekHealed > 4)  & (WeekHealed <= 8)) ~"<=8 Weeks",
  ((WeekHealed > 8)  & (WeekHealed <= 12)) ~"<=12 Weeks",
  (WeekHealed > 12)  ~">12 Weeks")
)

Plot_Data$HealingTime = factor(Plot_Data$HealingTime, levels=c("<=4 Weeks", "<=8 Weeks", "<=12 Weeks", ">12 Weeks"))
healingcolors = c( "yellow","goldenrod1", "darkorange1","firebrick4")

apeplot_healing <- ggtree::ggtree(Plot_Data, size=.15, layout = "rectangular") + geom_tiplab(size=2.5) + geom_point(size=2, aes(color=HealingTime)) + scale_color_manual(values=healingcolors) + geom_treescale(x=0,y=-2) + theme(legend.text=element_text(size=20), legend.title=element_text(size=20), legend.key.size=unit(1, "cm"), legend.position = c(.1, .5) )+ guides(color = guide_legend(override.aes = list(size = 4))) #+ ggtitle("Phylogeny Labeled by Subject_Timepoint")
apeplot_healing$data$label = apeplot_healing$data$SubjTime

apeplotBar = ggtree::ggtree(Plot_Data, size=1, layout = "rectangular") + geom_tiplab() 

apeplotCircle = ggtree::ggtree(apeplotBar$data, layout="circular") + geom_tiplab()
ggsave(apeplotCircle, file="CircularPreliminaryPlotSBDRC.pdf", width=25,height=25)

apeplotBar$data$label = apeplotBar$data$CC
apeplotBarForCollapse = apeplotBar  + geom_text2(aes(subset=!isTip, label=node)) +  geom_tiplab(size=2.5)

ggsave(apeplotBarForCollapse, file="NodesLabeledCutoffPreliminaryCCs.pdf", width=25, height=40)


NodesCollapse = c(298, 286, 253, 488, 491, 462, 474, 485, 484, 341, 433, 456,443)

# for getting descendants of each node 
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
barplotdata = barplotdata %>% mutate(HealingTime=case_when(
  (WeekHealed <= 4) ~"<=4 Weeks", 
  ((WeekHealed > 4)  & (WeekHealed <= 8)) ~"<=8 Weeks",
  ((WeekHealed > 8)  & (WeekHealed <= 12)) ~"<=12 Weeks",
  (WeekHealed > 12)  ~">12 Weeks")
)
barplotdata_weekstoheal = barplotdata[c("grouping","patient_id", "HealingTime")]
barplotdata = barplotdata %>% mutate(SubjTime = if_else(is.na(SubjTime), label, SubjTime))

barplotdata = barplotdata %>% mutate(Healedby12 = HealingTime %in% c("<=4 Weeks", "<=8 Weeks", "<=12 Weeks"))

barplotdata_HealedEver = barplotdata[c("grouping","patient_id", "HealedYorN")]
barplotdata_weekstoheal = barplotdata[c("grouping","patient_id", "HealingTime")]

barplotdata = barplotdata[c("grouping","patient_id", "Healedby12")]
barplotdata_melted = barplotdata %>% filter(grouping!=0) %>% filter(!is.na(patient_id))  %>% reshape2::melt(id.vars=c("grouping","patient_id" ))

barplotdata_melted_uniqueSubjects = barplotdata_melted %>% group_by(grouping, patient_id) %>% unique()


BarPlotNew  <- ggplot(barplotdata_melted_uniqueSubjects, aes(x=as.factor(grouping), fill=(value))) + geom_bar() + coord_flip() + scale_fill_manual(values=c("darkred", "#9FE2BF"))  + theme(axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size=40), legend.position="bottom", axis.line.x=element_line())+ scale_x_discrete(limits=rev(NodeOrder)) 
BarPlotPrint = BarPlotNew + theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background  = element_rect(fill = "transparent", colour = NA), axis.ticks.y = element_blank(), axis.text.y=element_blank())


ggsave(file="CollapsedSBDRCTree_labelednode.pdf", collapsedTree +geom_treescale(linesize = 2) , width=15,height=20)
ggsave(file="CollapsedSBDRCBarPlot_labelednode.pdf", BarPlotPrint, width=23,height=20)



empty =ggplot(barplotdata_melted_uniqueSubjects, aes(x=as.factor(grouping), fill=value)) + 
  geom_tile(aes(fill=value, y=)) + scale_fill_viridis_c() + 
  theme_minimal() + xlab(NULL) + ylab(NULL)

empty %>% insert_left(collapsedTree) %>% insert_right(BarPlotPrint, width=.5) 
cowplot::plot_grid(collapsedTree , BarPlotNew, ncol=2, rel_widths=c(1, .5))
# + geom_tiplab(align=TRUE)
unique(apeplotCircle$data$CC)

circdata = apeplotCircle$data %>% mutate(CCRefs =  case_when(label=="CC12" ~ "CC12", 
                                                  label=="CC133" ~ "CC133",
                                                  label=="CC15" ~ "CC15",
                                                  label=="CC1_MSSA476"~ "CC1",
                                                  label=="CC1_MW2" ~ "CC1",
                                                  label=="CC1-ST188" ~ "CC1",
                                                  label=="CC22_HO" ~ "CC22",
                                                  label=="CC30_MRSA252" ~ "CC30",
                                                  label=="CC398" ~ "CC398",
                                                  label=="CC45" ~ "CC45",
                                                  label=="CC59" ~ "CC59",
                                                  label=="CC5_Mu50"~ "CC5",
                                                  label=="CC5_N315"~ "CC5",
                                                  label=="CC72_CN1" ~ "CC72",
                                                  label=="CC7" ~ "CC7",
                                                  label=="CC80" ~ "CC80",
                                                  label=="CC8_NCTC8325" ~ "CC8",
                                                  label=="CC8_Newman" ~ "CC8",
                                                  label=="CC97_ATCC6538" ~ "CC97",
                                                  label=="CC97_Newbould305" ~ "CC97",
                                                  label=="SA_502A" ~ "CC5",
                                                  label=="SA_AR464" ~ "CC45",
                                                  label=="S_epidermidis" ~ "S_epidermidis",
                                                  label=="ST20"  ~ "CC20",
                                                  label=="USA100_AR465"~ "CC5",
                                                  label=="USA300_FPR3757" ~ "CC8",
                                                  label=="USA400_051" ~ "CC1",
                                                  TRUE ~ as.character(CC)))
randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#FFFF33",
                "#33A02C","#E6F5C9")
circdataNullLabel=circdata
circdataNullLabel$label=""
circdataNullLabel$TipLABEL=NULL
CircularPlotColors = ggtree(circdataNullLabel, layout="circular") + geom_fruit(aes(y=CCRefs,x=0,fill=CCRefs), geom=geom_tile,offset=5, size=.1)  + geom_tiplab(size=0)
CircularPlotColors = CircularPlotColors + scale_fill_manual(values=randpalette18)
ggsave(CircularPlotColors, file="ColorsCCplotSBDRC.pdf", width=26, height=25)


geom_fruit(data=dat2, geom=geom_tile,
           mapping=aes(y=ID, x=Sites, alpha=Abundance, fill=Sites),
           color = "grey50", offset = 0.04,size = 0.02)

ggsave(apeplot, file="VeryPreliminarySBDRC.pdf", width=25, height=25)

