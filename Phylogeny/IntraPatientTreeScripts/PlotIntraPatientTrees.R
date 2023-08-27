# Amy Campbell
# Intra-patient trees of all DFU containing at least 3 isolates from the same CC 
# From ClonalframeML-corrected newick files
library("dplyr")
library("ape")
library("ggtree")
library("dplyr")
library("ggplot2")
library("ggtreeExtra")
library("ggstar")
library("ggnewscale")
library("cowplot")
library("stringr")
library("RColorBrewer")
library("reshape2")
library("stats")
library("cluster")

#setwd('/Users/amycampbell/Documents/DataInputGithub/')

SakGenomes = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence.csv") %>% filter(Annotation == "staphylokinase")

SakGenomes = SakGenomes[15:ncol(SakGenomes)]
SakGenomes[SakGenomes==""] <- 0
SakGenomes[SakGenomes!=0] <- 1
SakGenomes = apply(SakGenomes, 2, as.numeric)
SakGenomes = colSums(SakGenomes)
GenomesWithSak=names(SakGenomes[SakGenomes>0])
GenomesWithSak = sapply(GenomesWithSak, function(x) str_replace(x, "DORN", "SA"))

StaphIsolateDORNs=read.csv("data/DFU_Staph_aureus_isolates.csv")
UpdatedPhenotypes=read.csv("~/Desktop/GriceLabGit/Staphyloxanthin/Data/InVitroData/Phenotypes_Data.csv")
StaphIsolateDORNs$IsolateID = paste0("SA", StaphIsolateDORNs$Doern.lab.bank.)
StaphIsolateDORNs$DORN = paste0("DORN", StaphIsolateDORNs$Doern.lab.bank.)
CCinfo = read.csv("data/Phylogeny2022Data/CCMapPlotting.csv")
CCinfo = CCinfo %>% left_join(StaphIsolateDORNs %>% select(DORN, patient_id), by="DORN")

phageContents = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/PresenceAbsenceHClustPhages.csv")
phageContents$IsolateID = sapply(phageContents$Genome,function(x) str_replace(x,"DORN", "SA"))
phageContents$Genome=NULL
phageContents$X = NULL

plasmidContents = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/Plasmid_Presence_Absence_UTD.csv")
plasmidContents$IsolateID = plasmidContents$X

plasmidContents$X=NULL
plasmidContents$X.1 = NULL


DistOutputPath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Distances/"


TotalPlas =rowSums(plasmidContents[1:(ncol(plasmidContents)-1)])
plasmidContents$IsolateID = sapply(plasmidContents$IsolateID, function(x) str_replace(x,"DORN", "SA"))
megadf =plasmidContents %>% left_join(phageContents, by="IsolateID")
allplasmid_phages = colnames(megadf)[colnames(megadf)!="IsolateID"]


HealingInfo = UpdatedPhenotypes %>% mutate(HealedBy12 = if_else(is.na(week_healed) | week_healed>12, "no", "yes")) %>% select(patient, HealedBy12) %>% unique()


# Set up color by week
#######################################################
WeekPalette =RColorBrewer::brewer.pal(11, "Spectral")

darkpurple =RColorBrewer::brewer.pal(11, "PuOr")[10]
WeekPalette = c("#8B0000", WeekPalette)

WeekPalette = c("black", WeekPalette)
WeekPalette = append(WeekPalette, darkpurple)
colormap=data.frame(colors=WeekPalette, week=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26))


# Min max normalize the phenotypes
###################################
Phenotypes = UpdatedPhenotypes %>% select(staphyloxanthin,staphylokinase, siderophore, biofilm, IsolateID, patient, visit, week_healed)
Phenotypes$staphyloxanthin = ((Phenotypes$staphyloxanthin  - min(Phenotypes$staphyloxanthin))/(max(Phenotypes$staphyloxanthin) - min(Phenotypes$staphyloxanthin)))

Phenotypes$staphylokinase = ((Phenotypes$staphylokinase  - min(Phenotypes$staphylokinase, na.rm=T))/(max(Phenotypes$staphylokinase, na.rm = T) - min(Phenotypes$staphylokinase,na.rm = T)))

Phenotypes$biofilm = ((Phenotypes$biofilm  - min(Phenotypes$biofilm))/(max(Phenotypes$biofilm) - min(Phenotypes$biofilm)))
Phenotypes$siderophore = ((Phenotypes$siderophore  - min(Phenotypes$siderophore))/(max(Phenotypes$siderophore) - min(Phenotypes$siderophore)))
Phenotypes$week_collected = Phenotypes$visit*2



CFTreeFolder="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/ClonalFrameTrees"
TreeFilesList = list.files(CFTreeFolder)

TreeFilesListPatientsRepresented = sapply(TreeFilesList, function(x) str_split(x, "_")[[1]][2])
PatientSummary = data.frame(patient=TreeFilesListPatientsRepresented)
HealingInfo$patient = sapply(HealingInfo$patient, as.character)

PatientSummary$patientCC = paste(PatientSummary$patient, PatientSummary$CC, sep="_")
CCinfo$patientCC = paste(CCinfo$patient_id, CCinfo$CCLabel, sep="_")
CCinfo %>% filter(patientCC %in% PatientSummary$patientCC )

write.csv(PatientSummary %>% left_join(HealingInfo, by="patient"), file="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/PlotSummaryIntraPatientTrees.csv")

PatientSummary$CC = sapply(TreeFilesList, function(x) str_split(x, "_")[[1]][3])

# HealingInfo %>% filter(HealedBy12 == "no")


TreeFile="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/ClonalFrameTrees/patient_176_CC5_clonalframeML.newick.labelled_tree.newick" #CC30_clonalframeML.newick.labelled_tree.newick"
TreeFileName = basename(TreeFile)
TreeString = str_replace(TreeFileName, "_clonalframeML.newick.labelled_tree.newick", "")
TreeString = str_replace(TreeString, ".newick", "")
cclist=c()


mytree = ggtree::read.tree(TreeFile)


# List of pairs of genomes flagged for follow-up
GenomePairDF = data.frame()

# List of trees flagged for follow-up
SmallDistPairList = c()


for(f in 1:length(TreeFilesList)){
  ####################################################################
  # Read in each clonalframeML-corrected intra-patient/within-CC tree
  ####################################################################
  TreeFilename = TreeFilesList[f]
  TreeFile=paste0(CFTreeFolder, "/",TreeFilename)
  TreeInput= ggtree::read.tree(TreeFile)
  TreeString = str_replace(TreeFilename, "_clonalframeML.newick.labelled_tree.newick", "")
  WhichCC=str_split(TreeString, "_")[[1]][3]
  cclist=append(cclist, WhichCC)
  TreeString = str_replace(TreeString, ".newick", "")
  
  
  genomenames = TreeInput$tip.label
  genomenames = genomenames[(grepl("CC", genomenames) | grepl("ST", genomenames))]
  RefGenome=genomenames[1]
  
  SortedLengths  = sort(TreeInput$edge.length)
  SecondLongest = SortedLengths[length(SortedLengths)-1]
  
  RootedTree = ape::root(TreeInput, RefGenome, resolve.root=T)
  
  TreePlot = ggtree::ggtree(RootedTree, size=.5, layout = "rectangular") + geom_tiplab(size=10)
  
  TreePlot$data[TreePlot$data$label %in% c(RefGenome), "x"] = SecondLongest
  
  TreePlot$data$IsolateID = TreePlot$data$label

  TreePlot$data = TreePlot$data %>% left_join(Phenotypes, by="IsolateID")
  
  # Get cophenetic distances and remove the reference genome which was just used for rooting
  CopheneneticDist = cophenetic(RootedTree)
  CopheneneticDist = CopheneneticDist[row.names(CopheneneticDist)[row.names(CopheneneticDist)!=RefGenome], ] 
  CopheneneticDist = CopheneneticDist[,colnames(CopheneneticDist)!=RefGenome]
  
  # make groups where minimum distance between any two points in different clusters is 1e-5
  PhyDistance = hclust(as.dist(CopheneneticDist), method="single")
  PhyClusterMembership = cutree(PhyDistance, h=1e-5)
  
  ClusterMem = data.frame(IsolateID = row.names(CopheneneticDist), Cluster=PhyClusterMembership)
  
  # Agglomerative clustering with average/mean linkage ; cluster if mean(distance) between clusters >.25
  for(cluster in unique(ClusterMem$Cluster)){
    IsolatesInCluster = (ClusterMem %>% filter(Cluster == cluster))$IsolateID
    
    if(length(IsolatesInCluster) > 1){
      ClusterIsolates =  TreePlot$data %>% filter(IsolateID %in% IsolatesInCluster)
      
      StaphyloxanthinDF = ClusterIsolates %>% select(IsolateID, staphyloxanthin) %>% filter(!is.na(staphyloxanthin))
      SiderophoreDF = ClusterIsolates%>% select(IsolateID, siderophore) %>% filter(!is.na(siderophore))
      BiofilmDF =ClusterIsolates %>% select(IsolateID, biofilm) %>% filter(!is.na(biofilm))
      StaphylokinaseDF = ClusterIsolates %>% select(IsolateID, staphylokinase) %>% filter(!is.na(staphylokinase))
      
      
      # Staphyloxanthin
      if(nrow(StaphyloxanthinDF)>1){
      
      StaphyloxanthinClusters = cutree(hclust(dist(StaphyloxanthinDF$staphyloxanthin), method="average"), h=.25)
      if(length(unique(StaphyloxanthinClusters))>1) {
        StaphyloxanthinDF$Cluster =StaphyloxanthinClusters
        write.csv(StaphyloxanthinDF %>% arrange(Cluster), paste0(DistOutputPath,"Comparisons_", TreeString,"_Staphyloxanthin_",cluster  , ".csv"  ))
      }
      }
      # Biofilm 
      
      if(nrow(BiofilmDF)>1){
      
      BiofilmClusters = cutree(hclust(dist(BiofilmDF$biofilm), method="average"), h=.25)
      if(length(unique(BiofilmClusters))>1) {
        BiofilmDF$Cluster =BiofilmClusters
        write.csv(BiofilmDF %>% arrange(Cluster), paste0(DistOutputPath,"Comparisons_", TreeString,"_Biofilm_",cluster , ".csv" ))
        
      }   }
      
      # Staphylokinase 
      if(nrow(StaphylokinaseDF) > 1){
        StaphylokinaseClusters = cutree(hclust(dist(StaphylokinaseDF$staphylokinase), method="average"), h=.25)
        if(length(unique(StaphylokinaseClusters))>1) {
          StaphylokinaseDF$Cluster = StaphylokinaseClusters
          write.csv(StaphylokinaseDF %>% arrange(Cluster), paste0(DistOutputPath,"Comparisons_", TreeString,"_Staphylokinase_",cluster  , ".csv"  ))
          
        }
        
        StaphylokinaseDF_JustSak = StaphylokinaseDF %>% filter(IsolateID %in% GenomesWithSak)
        if(nrow(StaphylokinaseDF_JustSak) > 1){
          StaphylokinaseClustersWithSak = cutree(hclust(dist(StaphylokinaseDF_JustSak$staphylokinase), method="average"), h=.25)
          
          if(length(unique(StaphylokinaseClustersWithSak))>1){
            StaphylokinaseDF_JustSak$Cluster = StaphylokinaseClustersWithSak
            write.csv(StaphylokinaseDF_JustSak %>% arrange(Cluster), paste0(DistOutputPath,"SakContainingComparisons_", TreeString,"_Staphylokinase_",cluster  , ".csv"  ))
            
          }
          }
      }
    
      
      # Siderophore 
      if(nrow(SiderophoreDF)>1){
        SiderophoreClusters = cutree(hclust(dist(SiderophoreDF$siderophore), method="average"), h=.25)
        if(length(unique(SiderophoreClusters))>1) {
          SiderophoreDF$Cluster = SiderophoreClusters
          write.csv(SiderophoreDF %>% arrange(Cluster), paste0(DistOutputPath,"Comparisons_", TreeString,"_Siderophore_",cluster  , ".csv"  ))
          
        }
        
        
      }
      
      
    }
    
   
  }
  


  
  colorpalettedf = colormap %>% filter(week %in%  TreePlot$data$week_collected)
  TreePlot = TreePlot+ geom_tippoint(shape=15,size=5, aes(color=factor(week_collected))) +theme(legend.position="none") + scale_color_manual(values=colorpalettedf$colors)
  genomesIncluded = TreeInput$tip.label
  genomesIncluded = genomesIncluded[genomesIncluded!=RefGenome]
  heightsave=6*(5/length(genomesIncluded))

  TreePlot$data$staphylokinase[is.na(TreePlot$data$staphylokinase)] <- as.numeric(1.0)
  TreePlot$data$staphyloxanthin[is.na(TreePlot$data$staphyloxanthin)] <- as.numeric(1.0)
  TreePlot$data$biofilm[is.na(TreePlot$data$biofilm)] <- as.numeric(1.0)
  TreePlot$data$siderophore[is.na(TreePlot$data$siderophore)] <- as.numeric(1.0)

  TreePlot = TreePlot + geom_treescale(width=1e-6) + new_scale_fill()+ geom_fruit(aes(x=(staphyloxanthin)),fill="#B8860B",axis.params=list(add.axis=TRUE, axis="x"), geom=geom_bar, stat="identity", orientation="y", offset=.2,pwidth = .2)+
    new_scale_fill()+ geom_fruit(aes(x=(staphylokinase)),fill="#2D9D92", geom=geom_bar, stat="identity", orientation="y",axis.params=list(add.axis=TRUE, axis="x"),offset=0,pwidth = .2)+
  new_scale_fill()+ geom_fruit(aes(x=(biofilm)),fill="#83B44B", geom=geom_bar, stat="identity", orientation="y",axis.params=list(add.axis=TRUE, axis="x"),offset=0,pwidth = .2)+
  new_scale_fill()+ geom_fruit(aes(x=(siderophore)),fill="#8B0000", geom=geom_bar, stat="identity", orientation="y",axis.params=list(add.axis=TRUE, axis="x"),offset=0,pwidth = .2)
  ggsave(TreePlot, file=paste0("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/ClonalFrameTreesAnnotated/", TreeString, "_tree.pdf"),height=heightsave, width=14.1)


  
  ###########################################
  ## Phage and plasmid presence by tree tip:
  ###########################################

  plasmidContentsPatient = plasmidContents %>% filter(IsolateID %in% genomesIncluded )
  plasmidContentsPatient = plasmidContentsPatient %>% select_if(function(x) any(x!=0))
  idorder=plasmidContentsPatient$IsolateID
  plasmidContentsPatient = plasmidContentsPatient %>% select_if(function(x) any(x==0)) # differentially present/absent plasmids only
  plasmidContentsPatient$IsolateID = idorder
  plasmidContentsPatient$type = "Plasmid"
  plasmidcontentsmelt= plasmidContentsPatient %>% melt(id.vars=c("IsolateID","type"))

  phageContentsPatient = phageContents %>% filter(IsolateID %in% genomesIncluded )
  phageContentsPatient = phageContentsPatient %>% select_if(function(x) any(x!=0))
  isolateorderphage = phageContentsPatient$IsolateID
  phageContentsPatient = phageContentsPatient %>% select_if(function(x) any(x==0))
  phageContentsPatient$IsolateID = isolateorderphage
  phageContentsPatient$type="Phage"
  phageContentsMelt = phageContentsPatient %>% melt(id.vars=c("IsolateID","type"))


  if( (ncol(phageContentsMelt)) > 2 & (ncol(plasmidcontentsmelt) > 2 )){
    phageContentsMelt$value[phageContentsMelt$value==1] <- 2
    HGT_Info = rbind(plasmidcontentsmelt, phageContentsMelt)
    starshapeobj=c(11,15)
    fillvals=c("white","darkblue","#024B30")
    colvals = c("white","white","#024B30")

  }
  if((ncol(phageContentsMelt) > 2) & (ncol(plasmidcontentsmelt)==2) ){
    phageContentsMelt$value[phageContentsMelt$value==1] <- 2
    HGT_Info=phageContentsMelt
    starshapeobj=c(11)
    fillvals=c("white","#024B30")
    colvals = c("white","#024B30")
  }
  if((ncol(phageContentsMelt)==2) & (ncol(plasmidcontentsmelt)>2 )){
    HGT_Info=plasmidcontentsmelt
    fillvals=c("white","darkblue")
    colvals = c("white", "white")
    starshapeobj=c(15)

  }
  if( (ncol(phageContentsMelt)==2) & (ncol(plasmidcontentsmelt)==2)){
    HGT_Info=data.frame()
    print(paste0("No variably present/absent phages or plasmids in ",TreeString ))
  }

  if(nrow(HGT_Info) > 0){

    yorder = (TreePlot$data %>% arrange(y) %>% filter(grepl(label, pattern="SA") & (label!=RefGenome)))$label


    widthobj=10*(length(unique(HGT_Info$variable))/ 9)
    HGTplot = ggplot(HGT_Info, aes(y=IsolateID, x=variable, color=factor(value), fill=factor(value), starshape=factor(type)))+ geom_star(size=1,size=1) +scale_color_manual(values=colvals)+ scale_starshape_manual(values=starshapeobj)+
     scale_fill_manual(values=fillvals) + theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), axis.text.y=element_text(size=2),
                                                legend.position="none", axis.line.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                                                panel.background = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())+scale_x_discrete(position = "top") #+ coord_fixed(ratio=1)

    ThingOrder = sapply(unique(HGTplot$data$variable),as.character)
    HGTplot$data$IsolateID = factor(HGTplot$data$IsolateID, levels=c(yorder, setdiff(megadf$IsolateID, yorder)))
    HGTplot$data$variable = factor(HGTplot$data$variable, levels=c(ThingOrder, setdiff(allplasmid_phages, ThingOrder)))


    HGTplot = HGTplot + scale_y_discrete(limits= c(yorder, setdiff(megadf$IsolateID, yorder)), drop=FALSE) + scale_x_discrete(limits= c(ThingOrder, setdiff(allplasmid_phages, ThingOrder)), drop=FALSE)

    ggsave(HGTplot, file=paste0("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/ClonalFrameTreesAnnotated/", TreeString, "_HGTplot.pdf"),height=12.5,width=7)

  }
}  



