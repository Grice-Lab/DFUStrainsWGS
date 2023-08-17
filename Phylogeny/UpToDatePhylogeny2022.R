# amy campbell 
# Phylogeny of DORN isolates + their references (247 genomes total)
# Updated 05-2022
setwd('/Users/amycampbell/Documents/DataInputGithub/')
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
library("ggstar")
library("ggnewscale")
library("cowplot")
TreeFilePath = "data/Phylogeny2022Data/RAxML_bestTree.RaxMLTree2022.newick"
clonalframetreepath="data/Phylogeny2022Data/clonalframe_tree.newick.labelled_tree.newick"
CCmapping="data/Phylogeny2022Data/CCMapPlotting.csv"

StaphIsolateDORNs=read.csv("data/DFU_Staph_aureus_isolates.csv")
UpdatedPhenotypes=read.csv("data/staphyloxanthin_paper_data.csv")
PhenotypesUpdatedSak= read.csv("data/staphyloxanthin_paper_UpdatedSakClassifications.csv")
PhenotypesUpdatedSak = PhenotypesUpdatedSak %>% select(DORN, SakUpdated)
UpdatedPhenotypes$DORN=paste0("DORN", UpdatedPhenotypes$DORN)

StaphIsolateDORNs$DORN=paste0("DORN", StaphIsolateDORNs$Doern.lab.bank.)
CCgroupings=read.csv(CCmapping)

colorscheme_CCs = read.csv("~/Documents/DataInputGithub/data/Phylogeny2022Data/CC_Colorscheme.csv")


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


# Adding min-max-normalized phenotypes to each tree 

Phenotypes = UpdatedPhenotypes %>% select(staphyloxanthin,staphylokinase, siderophore, biofilm, DORN)
Phenotypes$staphyloxanthin = ((Phenotypes$staphyloxanthin  - min(Phenotypes$staphyloxanthin))/(max(Phenotypes$staphyloxanthin) - min(Phenotypes$staphyloxanthin)))

Phenotypes$staphylokinase = ((Phenotypes$staphylokinase  - min(Phenotypes$staphylokinase, na.rm=T))/(max(Phenotypes$staphylokinase, na.rm = T) - min(Phenotypes$staphylokinase,na.rm = T)))

Phenotypes$biofilm = ((Phenotypes$biofilm  - min(Phenotypes$biofilm))/(max(Phenotypes$biofilm) - min(Phenotypes$biofilm)))
Phenotypes$siderophore = ((Phenotypes$siderophore  - min(Phenotypes$siderophore))/(max(Phenotypes$siderophore) - min(Phenotypes$siderophore)))


# CC vs. each phenotype (kruskall-wallis)
#########################################

# Biofilm
############
PhenotypesCC = PhenotypesCC %>% filter(!is.na(CCLabel))
meansbiofilm= PhenotypesCC %>% select(CCLabel, biofilm) %>% group_by(CCLabel) %>% summarize(CCLabel=CCLabel, biofilm=mean(biofilm)) %>% arrange(biofilm) 
OrderBiofilm = unique(meansbiofilm$CCLabel)
kruskal.test(biofilm ~ CCLabel, data = PhenotypesCC)
pairwise.wilcox.test(PhenotypesCC$biofilm, PhenotypesCC$CCLabel,
                     p.adjust.method = "BH")
BiofilmPlot = ggplot(PhenotypesCC, aes(x=CCLabel, y=biofilm)) + geom_boxplot(fill="#83B44B") + xlab("Clonal Complex") + ylab("Normalized Biofilm Formation") + theme_classic() + ggpubr::stat_compare_means(method="kruskal.test", label.y=1.2, label.x=2.5) + ggtitle("Biofilm") + theme(plot.title=element_text(hjust=.5, face="bold", size=14))
BiofilmPlot$data$CCLabel = factor(BiofilmPlot$data$CCLabel, levels=OrderBiofilm)


# Staphyloxanthin
################
PhenotypesCC = PhenotypesCC %>% filter(!is.na(CCLabel))
kruskal.test(staphyloxanthin ~ CCLabel, data = PhenotypesCC)
pairwise.wilcox.test(PhenotypesCC$staphyloxanthin, PhenotypesCC$CCLabel,
                     p.adjust.method = "BH")
staphyloxanthinPlot = ggplot(PhenotypesCC, aes(x=CCLabel, y=staphyloxanthin)) + geom_boxplot(fill="#B8860B") + xlab("Clonal Complex") + ylab("Normalized Staphyloxanthin Production") + theme_classic() + ggpubr::stat_compare_means(method="kruskal.test", label.y=1.2, label.x=2.5) + ggtitle("Staphyloxanthin") + theme(plot.title=element_text(hjust=.5, face="bold", size=14))
staphyloxanthinPlot$data$CCLabel = factor(staphyloxanthinPlot$data$CCLabel, levels=OrderBiofilm)





# staphylokinase
################
PhenotypesCC = PhenotypesCC %>% filter(!is.na(CCLabel))
kruskal.test(staphylokinase ~ CCLabel, data = PhenotypesCC)
pairwise.wilcox.test(PhenotypesCC$staphylokinase, PhenotypesCC$CCLabel,
                     p.adjust.method = "BH")
staphylokinasePlot = ggplot(PhenotypesCC, aes(x=CCLabel, y=staphylokinase)) + geom_boxplot(fill="#2D9D92") + xlab("Clonal Complex") + ylab("Normalized Staphylokinase Production") + theme_classic() + ggpubr::stat_compare_means(method="kruskal.test", label.y=1.2, label.x=2.5) + ggtitle("Staphylokinase") + theme(plot.title=element_text(hjust=.5, face="bold", size=14))
staphylokinasePlot$data$CCLabel = factor(staphylokinasePlot$data$CCLabel, levels=OrderBiofilm)



# staphylokinase No sak
#######################
PhenotypesCCSak = PhenotypesCC %>% filter(!is.na(CCLabel)) %>% left_join((PhenotypesUpdatedSak %>% select(DORN, SakUpdated)),by="DORN" ) %>% filter(SakUpdated=="yes")
kruskal.test(staphylokinase ~ CCLabel, data = PhenotypesCCSak)
pairwise.wilcox.test(PhenotypesCCSak$staphylokinase, PhenotypesCCSak$CCLabel,
                     p.adjust.method = "BH")
staphylokinasePlotNoSak = ggplot(PhenotypesCCSak, aes(x=CCLabel, y=staphylokinase)) + geom_boxplot(fill="#2D9D92") + xlab("Clonal Complex") + ylab("Normalized Staphylokinase Production") + theme_classic() + ggpubr::stat_compare_means(method="kruskal.test", label.y=1.2, label.x=2.5) + ggtitle("Staphylokinase(sak-positive Only)") + theme(plot.title=element_text(hjust=.5, face="bold", size=14))
staphylokinasePlotNoSak$data$CCLabel = factor(staphylokinasePlotNoSak$data$CCLabel, levels=OrderBiofilm)


# Siderophore production
#######################
PhenotypesCC = PhenotypesCC %>% filter(!is.na(CCLabel))
kruskal.test(siderophore ~ CCLabel, data = PhenotypesCC)
pairwise.wilcox.test(PhenotypesCC$siderophore, PhenotypesCC$CCLabel,
                     p.adjust.method = "BH")
siderophorePlot = ggplot(PhenotypesCC, aes(x=CCLabel, y=siderophore)) + geom_boxplot(fill="#8B0000")  + geom_jitter(height=0) + xlab("Clonal Complex") + ylab("Normalized Siderophore Production") + theme_classic() + ggpubr::stat_compare_means(method="kruskal.test", label.y=1.2, label.x=2.5) + ggtitle("Siderophore") + theme(plot.title=element_text(hjust=.5, face="bold", size=14))
siderophorePlot$data$CCLabel = factor(siderophorePlot$data$CCLabel, levels=OrderBiofilm)


gridExtra::grid.arrange(staphyloxanthinPlot, BiofilmPlot, staphylokinasePlot, siderophorePlot)
gridExtra::grid.arrange(staphyloxanthinPlot, BiofilmPlot, staphylokinasePlotNoSak, siderophorePlot)


##############################################
# Line plots for phenotypes vs. CC (Figure 2A)
##############################################
# Phenotypes by CC 



Phenotypes_With_CC = UpdatedPhenotypes %>% left_join(CCgroupings, by="DORN") %>% filter(!is.na(CCLabel))
MeansByCC = Phenotypes_With_CC %>% group_by(CCLabel) %>% summarize(meanSiderophore=mean(siderophore, na.rm=T), meanStaphylokinase = mean(staphylokinase, na.rm=T), meanStaphyloxanthin=mean(staphyloxanthin, na.rm=T), meanBiofilm=mean(biofilm, na.rm=T))
SDsByCC = Phenotypes_With_CC %>% group_by(CCLabel) %>% summarize(sdSiderophore=sd(siderophore, na.rm=T), sdStaphylokinase = sd(staphylokinase, na.rm=T), sdStaphyloxanthin=sd(staphyloxanthin, na.rm=T), sdBiofilm=sd(biofilm, na.rm=T))
SDsByCC[is.na(SDsByCC)] <- 0


# staphyloxanthin line plot
MeansByCC$LowerXanthin = MeansByCC$meanStaphyloxanthin - SDsByCC$sdStaphyloxanthin
MeansByCC$UpperXanthin = MeansByCC$meanStaphyloxanthin + SDsByCC$sdStaphyloxanthin
MeansByCC = MeansByCC %>% left_join(colorscheme_CCs,by="CCLabel")
CCxanthinLinePlot = ggplot(MeansByCC,aes( x=factor(CCLabel), y=meanStaphyloxanthin))+ geom_point(color="#B8860B") + geom_errorbar(color="#B8860B",width=.2,aes(ymax=UpperXanthin, ymin=LowerXanthin)) 
CCxanthinLinePlot$data$CCLabel = factor(CCxanthinLinePlot$data$CCLabel, levels=rev(unique(CCxanthinLinePlot$data$CCLabel)))
CCxanthinLinePlot = CCxanthinLinePlot + theme_classic() + coord_flip() + theme(axis.text.y=element_text(colour=rev(MeansByCC$hexval)))
ggsave(CCxanthinLinePlot,file="~/Documents/Saureus_Genomics_Paper/Figure2Figs/STX_lineplot.pdf",height=6,width=2)
MeansByCC$meanStaphyloxanthin

# biofilm line plot 
MeansByCC$LowerBiofilm = MeansByCC$meanBiofilm - SDsByCC$sdBiofilm
MeansByCC$UpperBiofilm = MeansByCC$meanBiofilm + SDsByCC$sdBiofilm
CCbiofilmLinePlot = ggplot(MeansByCC,aes( x=factor(CCLabel), y=meanBiofilm))+ geom_point(color="#83B44B") + geom_errorbar(color="#83B44B",width=.2,aes(ymax=UpperBiofilm, ymin=LowerBiofilm)) 
CCbiofilmLinePlot$data$CCLabel = factor(CCbiofilmLinePlot$data$CCLabel, levels=rev(unique(CCbiofilmLinePlot$data$CCLabel)))
CCbiofilmLinePlot = CCbiofilmLinePlot + theme_classic() + coord_flip() + theme(axis.text.y=element_text(colour=rev(MeansByCC$hexval)))
ggsave(CCbiofilmLinePlot,file="~/Documents/Saureus_Genomics_Paper/Figure2Figs/Biofilm_lineplot.pdf",height=6,width=2)

xanthin_line = ggplot(MeansByCC) + aes(x=CCLabel,y)
 
# biofilm line plot 
MeansByCC$LowerSiderophore = MeansByCC$meanSiderophore - SDsByCC$sdSiderophore
MeansByCC$UpperSiderophore = MeansByCC$meanSiderophore + SDsByCC$sdSiderophore
CCsiderophoreLinePlot = ggplot(MeansByCC,aes( x=factor(CCLabel), y=meanSiderophore))+ geom_point(color="#8B0000") + geom_errorbar(color="#8B0000",width=.2,aes(ymax=UpperSiderophore, ymin=LowerSiderophore)) 
CCsiderophoreLinePlot$data$CCLabel = factor(CCsiderophoreLinePlot$data$CCLabel, levels=rev(unique(CCsiderophoreLinePlot$data$CCLabel)))
CCsiderophoreLinePlot = CCsiderophoreLinePlot + theme_classic() + coord_flip() + theme(axis.text.y=element_text(colour=rev(MeansByCC$hexval)))
ggsave(CCsiderophoreLinePlot,file="~/Documents/Saureus_Genomics_Paper/Figure2Figs/Siderophore_lineplot.pdf",height=6,width=2)

# biofilm line plot 
MeansByCC$LowerStaphylokinase = MeansByCC$meanStaphylokinase - SDsByCC$sdStaphylokinase
MeansByCC$UpperStaphylokinase = MeansByCC$meanStaphylokinase + SDsByCC$sdStaphylokinase

CCstaphylokinaseLinePlot = ggplot(MeansByCC,aes( x=factor(CCLabel), y=meanStaphylokinase))+ geom_point(color="#2D9D92") + geom_errorbar(color="#2D9D92",width=.2,aes(ymax=UpperStaphylokinase, ymin=LowerStaphylokinase)) 
CCstaphylokinaseLinePlot$data$CCLabel = factor(CCstaphylokinaseLinePlot$data$CCLabel, levels=rev(unique(CCstaphylokinaseLinePlot$data$CCLabel)))
CCstaphylokinaseLinePlot = CCstaphylokinaseLinePlot + theme_classic() + coord_flip() + theme(axis.text.y=element_text(colour=rev(MeansByCC$hexval)))
ggsave(CCstaphylokinaseLinePlot,file="~/Documents/Saureus_Genomics_Paper/Figure2Figs/Staphylokinase_lineplot.pdf",height=6,width=2)






apeplotcircPhenotypes<- ggtree::ggtree(apeRootedCF, size=.25, layout = "circular", branch.length = "none") #+ geom_tiplab(size=1.7)

apeplotcircPhenotypes$data$DORN = apeplotcirc$data$label
#apeplotcircPhenotypes$data$label = NA
apeplotcircPhenotypes$data = apeplotcircPhenotypes$data %>% left_join(CCgroupings, by="DORN")
apeplotcircPhenotypes$data  = apeplotcircPhenotypes$data %>% mutate(CCRefs =  case_when(label=="CC1_MSSA476" ~"CC1",
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




apeplotcircPhenotypes$data = apeplotcircPhenotypes$data %>% left_join(Phenotypes, by="DORN")
apeplotcircPhenotypes$data = apeplotcircPhenotypes$data %>% left_join(StaphIsolateDORNs, by="DORN")
apeplotcircPhenotypes$data$label=apeplotcircPhenotypes$data$patient_id

# Plot with bars for staphyloxanthin production & CC 
#####################################################
WithXanthin = apeplotcircPhenotypes + geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=4, pwidth=6, color=NA, offset=.08) +
  scale_fill_manual(values=randpalette18)+  new_scale_fill()+
  geom_fruit(aes(x=(staphyloxanthin)),fill="#B8860B", geom=geom_bar, stat="identity", orientation="y",offset=.1, pwidth=10)
commonlegend = get_legend(WithXanthin)

ggsave(WithXanthin, file="data/Phylogeny2022Data/TestXanthinTree.pdf", width=11,height=11)
WithXanthin = WithXanthin + theme(legend.position="None")

# Plot with bars for staphylokinase production & CC 
#####################################################
WithKinase = apeplotcircPhenotypes + geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=4, pwidth=6, color=NA) +
  scale_fill_manual(values=randpalette18)+  new_scale_fill()+
  geom_fruit(aes(x=(staphylokinase)),fill="#2D9D92", geom=geom_bar, stat="identity", orientation="y",offset=.1, pwidth=10)

ggsave(WithKinase, file="data/Phylogeny2022Data/TestKinaseTree.pdf", width=11,height=11)


WithKinaseAndGene = apeplotcircPhenotypes + geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=4, pwidth=6, color=NA) +
  scale_fill_manual(values=randpalette18)+  new_scale_fill()+
  geom_fruit(aes(x=(staphylokinase)),fill="#2D9D92", geom=geom_bar, stat="identity", orientation="y",offset=.1, pwidth=10)

WithKinaseAndGene$data = WithKinaseAndGene$data %>% left_join(PhenotypesUpdatedSak, by="DORN")
WithKinaseAndGene$data$SakUpdated = if_else(is.na(WithKinaseAndGene$data$SakUpdated), "Reference", WithKinaseAndGene$data$SakUpdated)
WithKinaseAndGene = WithKinaseAndGene +geom_tippoint(aes(color=factor(SakUpdated)), shape=16,position=position_nudge(x=1)) + scale_color_manual(values=(c("transparent", "transparent","black")))

ggsave(WithKinaseAndGene, file="data/Phylogeny2022Data/TestKinaseTreeWithGene.pdf", width=11,height=11)


WithKinase = WithKinase + theme(legend.position="None")



# Plot with bars for siderophore production & CC
################################################
WithSiderophore = apeplotcircPhenotypes + geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=4, pwidth=6, color=NA, offset=.08) +
  scale_fill_manual(values=randpalette18)+  new_scale_fill()+
  geom_fruit(aes(x=(siderophore)),fill="#8B0000", geom=geom_bar, stat="identity", orientation="y",offset=.1, pwidth=10)
ggsave(WithSiderophore, file="data/Phylogeny2022Data/TestSiderophoreTree.pdf", width=11,height=11)

WithSiderophore= WithSiderophore + theme(legend.position="None")


# Plot with bars for siderophore production & CC
################################################
WithBiofilm = apeplotcircPhenotypes + geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=4, pwidth=6, color=NA, offset=.08) +
  scale_fill_manual(values=randpalette18)+  new_scale_fill()+
  geom_fruit(aes(x=(biofilm)),fill="#83B44B", geom=geom_bar, stat="identity", orientation="y",offset=.1, pwidth=10)

ggsave(WithBiofilm, file="data/Phylogeny2022Data/TestBiofilmTree.pdf", width=11,height=11)
WithBiofilm = WithBiofilm + theme(legend.position="None")

# WithBiofilm
##################################

grid1  = plot_grid(WithXanthin, WithKinase, WithSiderophore,WithBiofilm )

##############################################################
# Full circular plot where each isolate is labeled by healing
##############################################################
circdata_healing = ggtree(circdata, layout="circular", size=.5) + geom_tiplab(size=0) +
  geom_tippoint(aes(color = HealedBy12),  size=2) + scale_color_manual(values=c("darkred", "#9FE2BF"),na.value ="#00000000" ) +
  new_scale_fill()+ geom_fruit(geom=geom_star, aes(y=label,fill=CCRefs), starshape=13,size=8,color=NA, offset=.05) + scale_fill_manual(values=randpalette18)

circdata_healing


ggsave(circdata_healing, file="CircularPlotCF_HealingAlone.pdf", width=10, height=10)
ggsave(circdata_healing, file="CircularPlotCF_Healing.pdf", width=12, height=12)


#############################################################################
# Make simplified version of the full rectangular tree with healing outcomes
##############################################################################
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

