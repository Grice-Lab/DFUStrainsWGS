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

setwd('/Users/amycampbell/Documents/DataInputGithub/')

StaphIsolateDORNs=read.csv("data/DFU_Staph_aureus_isolates.csv")
UpdatedPhenotypes=read.csv("~/Desktop/GriceLabGit/Staphyloxanthin/Data/InVitroData/Phenotypes_Data.csv")
StaphIsolateDORNs$IsolateID = paste0("SA", StaphIsolateDORNs$Doern.lab.bank.)

phageContents = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/PhagePresenceAbsence.csv")
phageContents$IsolateID = sapply(phageContents$X,function(x) str_replace(x,"DORN", "SA"))
phageContents$X = NULL

plasmidContents = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/Filtered_Plasmid_Presence_Absence.csv")
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



write.csv(PatientSummary %>% left_join(HealingInfo, by="patient"), file="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/PlotSummaryIntraPatientTrees.csv")

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
  print(TreeString)
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
  
  CopheneneticDist = cophenetic(RootedTree)
  CopheneneticDist = CopheneneticDist[row.names(CopheneneticDist)[row.names(CopheneneticDist)!=RefGenome], ] 
  CopheneneticDist = CopheneneticDist[,colnames(CopheneneticDist)!=RefGenome]
  
  DistanceDF = reshape2::melt(CopheneneticDist)
  colnames(DistanceDF) = c("Genome1", "Genome2", "CopheneticDist")
  


  StaphyloxanthinDF = TreePlot$data %>% select(IsolateID, staphyloxanthin) %>% filter(!is.na(staphyloxanthin))
  SiderophoreDF = TreePlot$data %>% select(IsolateID, siderophore) %>% filter(!is.na(siderophore))
  BiofilmDF = TreePlot$data %>% select(IsolateID, biofilm) %>% filter(!is.na(biofilm))
  StaphylokinaseDF = TreePlot$data %>% select(IsolateID, staphylokinase) %>% filter(!is.na(staphylokinase))

  
  StaphylokinaseDist = data.frame(as.matrix(dist(StaphylokinaseDF$staphylokinase,diag=T, upper=T)))
  colnames(StaphylokinaseDist)= StaphylokinaseDF$IsolateID
  StaphylokinaseDist$IsolateID = StaphylokinaseDF$IsolateID
  MeltedStaphylokinase = reshape2::melt(StaphylokinaseDist)
  colnames(MeltedStaphylokinase) = c("Genome1", "Genome2", "StaphylokinaseDistance")
  

  
  StaphyloxanthinDist = data.frame(as.matrix(dist(StaphyloxanthinDF$staphyloxanthin,diag=T, upper=T)))
  colnames(StaphyloxanthinDist)= StaphyloxanthinDF$IsolateID
  StaphyloxanthinDist$IsolateID = StaphyloxanthinDF$IsolateID
  MeltedStaphyloxanthin = reshape2::melt(StaphyloxanthinDist)
  colnames(MeltedStaphyloxanthin) = c("Genome1", "Genome2", "StaphyloxanthinDistance")
  

  
  
  
  SiderophoreDist = data.frame(as.matrix(dist(SiderophoreDF$siderophore,diag=T, upper=T)))
  colnames(SiderophoreDist)= SiderophoreDF$IsolateID
  SiderophoreDist$IsolateID = SiderophoreDF$IsolateID
  MeltedSiderophore = reshape2::melt(SiderophoreDist)
  colnames(MeltedSiderophore) = c("Genome1", "Genome2", "SiderophoreDistance")
  
  
  BiofilmDist = data.frame(as.matrix(dist(BiofilmDF$biofilm,diag=T, upper=T)))
  colnames(BiofilmDist)= BiofilmDF$IsolateID
  BiofilmDist$IsolateID = BiofilmDF$IsolateID
  MeltedBiofilm = reshape2::melt(BiofilmDist)
  colnames(MeltedBiofilm) = c("Genome1", "Genome2", "BiofilmDistance")
  
  
  AllDistances = DistanceDF %>% left_join(MeltedStaphylokinase,by=c("Genome1", "Genome2")) %>%
    left_join(MeltedStaphyloxanthin, by=c("Genome1", "Genome2")) %>%
    left_join(MeltedSiderophore, by=c("Genome1", "Genome2")) %>%
    left_join(MeltedBiofilm, by=c("Genome1", "Genome2")) 
  
  
  
  
  NeedToCompareStaphyloxanthin = AllDistances %>% filter((CopheneticDist < 1e-05) & (StaphyloxanthinDistance > .2))
  NeedToCompareStaphylokinase = AllDistances %>% filter((CopheneticDist < 1e-05) & (StaphylokinaseDistance > .2))
  NeedToCompareSiderophore = AllDistances %>% filter((CopheneticDist < 1e-05) & (SiderophoreDistance > .2))
  NeedToCompareBiofilm = AllDistances %>% filter((CopheneticDist < 1e-05) & (BiofilmDistance > .2))
  
  # Biofilm comparisons
  #####################
  print("Biofilm")
  
  ComparisonsList = c()
  if(nrow(NeedToCompareBiofilm) > 0 ){
    genomes_to_compare = unique((NeedToCompareBiofilm %>% arrange(-BiofilmDistance))$Genome1)
    
    i=1
    while(length(genomes_to_compare)>0){
      currentgenome = genomes_to_compare[1]
      DistSubset = AllDistances %>% filter(Genome1==currentgenome)
      
      A_Phenotype_Genomes = (DistSubset %>% filter( (BiofilmDistance <= .2) & (CopheneticDist < 1e-5)))$Genome2
      B_Phenotype_Genomes = (DistSubset %>% filter( (BiofilmDistance > .2) & (CopheneticDist < 1e-5)))$Genome2
      
   
      # check that genomes in two different groups at least have .1 distance 
      CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
      CheckDiffsGenomes = paste(CheckDiffs$Genome1, CheckDiffs$Genome2, sep="_")
      
      ChooseToExclude =  (CheckDiffs[which(CheckDiffs$BiofilmDistance <= .2),])
      while(nrow(ChooseToExclude)>0){
        
   
        # if there's a kind of 'edge case' where there's some match of A vs. B genome without much phenotype distance (not >.2)
        # Remove the one with the maximum cophenetic distance from other genomes that's the greatest
        ListPotentialExcludeA = sapply(ChooseToExclude$Genome1, as.character)
        ListPotentialExcludeB = sapply(ChooseToExclude$Genome2, as.character)
        
        
        if((length(ListPotentialExcludeA) + length(ListPotentialExcludeB)) > 0){
          PotentialExcludes = c(ListPotentialExcludeA, ListPotentialExcludeB)
          if(length(B_Phenotype_Genomes) ==1){
            PotentialExcludes = ListPotentialExcludeA
          }
          if(length(A_Phenotype_Genomes)==1){
            PotentialExcludes = ListPotentialExcludeB
          }
            maxdist = 0
            for(p in PotentialExcludes){
              current = max( ( AllDistances %>% filter(Genome1==p))$CopheneticDist  ) 
              if(current > maxdist){
                remove=p
                maxdist=current
              }
                
            }
          }
        B_Phenotype_Genomes = B_Phenotype_Genomes[B_Phenotype_Genomes!=remove]
        A_Phenotype_Genomes = A_Phenotype_Genomes[A_Phenotype_Genomes!=remove]
        CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
        CheckDiffsGenomes = paste(CheckDiffs$Genome1, CheckDiffs$Genome2, sep="_")
        
        ChooseToExclude =  (CheckDiffs[which(CheckDiffs$BiofilmDistance < .2),])
        
      }
      A_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% A_Phenotype_Genomes))
      B_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% B_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes))
      
      ACheck= ( (max(A_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(A_Phenotype_GenomeComparisons$BiofilmDistance) <=2 ) )
      BCheck= ( (max(B_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(B_Phenotype_GenomeComparisons$BiofilmDistance) <=2 ) )
      
      if(length(ComparisonsList)>0){
        LengthsIntersect = sapply(1:length(ComparisonsList),function(x) length(intersect(ComparisonsList[[x]], CheckDiffsGenomes)))
      }else{
        LengthsIntersect=c()
      }
      
      if(ACheck==T & BCheck==T & length(B_Phenotype_Genomes)>0 & length(A_Phenotype_Genomes)>0 & !any(LengthsIntersect>0)){
        PhenotypesA = TreePlot$data %>% select(IsolateID, biofilm) %>% filter(IsolateID %in% A_Phenotype_Genomes) 
        PhenotypesA$Phenotype="A"
        
        PhenotypesB = TreePlot$data %>% select(IsolateID, biofilm) %>% filter(IsolateID %in% B_Phenotype_Genomes) 
        PhenotypesB$Phenotype="B"
        
        result=rbind(PhenotypesA, PhenotypesB)
        write.csv(result,paste0(DistOutputPath, TreeString, "_BiofilmComparisons",toString(i), ".csv"))
        ComparisonsList[[i]] = CheckDiffsGenomes
        
        genomes_to_compare = genomes_to_compare[!(genomes_to_compare %in% result$IsolateID) ]
        i=i+1
        
      }else{
        genomes_to_compare = genomes_to_compare[genomes_to_compare!=currentgenome]
      }
    }}
  
      
  
  
  # Siderophore comparisons
  #####################
  print("Siderophore")
  
  ComparisonsList = c()
  i=1
  
  if(nrow(NeedToCompareSiderophore) > 0 ){
    
    genomes_to_compare = unique((NeedToCompareSiderophore %>% arrange(-SiderophoreDistance))$Genome1)
    
    while(length(genomes_to_compare)>0){
      currentgenome = genomes_to_compare[1]
      DistSubset = AllDistances %>% filter(Genome1==currentgenome)
      
      A_Phenotype_Genomes = (DistSubset %>% filter( (SiderophoreDistance <= .2) & (CopheneticDist < 1e-5)))$Genome2
      B_Phenotype_Genomes = (DistSubset %>% filter( (SiderophoreDistance > .2) & (CopheneticDist < 1e-5)))$Genome2
      

      # check that genomes in two different groups at least have .1 distance 
      CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
      CheckDiffsGenomes = paste(CheckDiffs$Genome1, CheckDiffs$Genome2, sep="_")
      
      ChooseToExclude =  (CheckDiffs[which(CheckDiffs$SiderophoreDistance <= .2),])
      while(nrow(ChooseToExclude)>0){
        
        # if there's a kind of 'edge case' where there's some match of A vs. B genome without much phenotype distance (not >.2)
        # Remove the one with the maximum cophenetic distance from other genomes that's the greatest
        ListPotentialExcludeA = sapply(ChooseToExclude$Genome1, as.character)
        ListPotentialExcludeB = sapply(ChooseToExclude$Genome2, as.character)
        
        
        if((length(ListPotentialExcludeA) + length(ListPotentialExcludeB)) > 0 ){
          PotentialExcludes = c(ListPotentialExcludeA, ListPotentialExcludeB)
          if(length(B_Phenotype_Genomes) ==1){
            PotentialExcludes = ListPotentialExcludeA
          }
          if(length(A_Phenotype_Genomes)==1){
            PotentialExcludes = ListPotentialExcludeB
          }
          maxdist = 0
          for(p in PotentialExcludes){
            current = max( ( AllDistances %>% filter(Genome1==p))$CopheneticDist  ) 
            if(current > maxdist){
              remove=p
              maxdist=current
            }
            
          }
        }
        B_Phenotype_Genomes = B_Phenotype_Genomes[B_Phenotype_Genomes!=remove]
        A_Phenotype_Genomes = A_Phenotype_Genomes[A_Phenotype_Genomes!=remove]
        CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
        CheckDiffsGenomes = paste(CheckDiffs$Genome1, CheckDiffs$Genome2, sep="_")
        
        ChooseToExclude =  (CheckDiffs[which(CheckDiffs$SiderophoreDistance < .2),])
        
      }
      
      if(length(ComparisonsList)>0){
        LengthsIntersect = sapply(1:length(ComparisonsList),function(x) length(intersect(ComparisonsList[[x]], CheckDiffsGenomes)))
      }else{
        LengthsIntersect=c()
      }
      
      A_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% A_Phenotype_Genomes))
      B_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% B_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes))
      
      ACheck= ( (max(A_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(A_Phenotype_GenomeComparisons$SiderophoreDistance) <=2 ) )
      BCheck= ( (max(B_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(B_Phenotype_GenomeComparisons$SiderophoreDistance) <=2 ) )
      
      if(ACheck==T & BCheck==T & length(B_Phenotype_Genomes)>0 & length(A_Phenotype_Genomes)>0 & !any(LengthsIntersect>0)){
        PhenotypesA = TreePlot$data %>% select(IsolateID, siderophore) %>% filter(IsolateID %in% A_Phenotype_Genomes) 
        PhenotypesA$Phenotype="A"
        
        PhenotypesB = TreePlot$data %>% select(IsolateID, siderophore) %>% filter(IsolateID %in% B_Phenotype_Genomes) 
        PhenotypesB$Phenotype="B"
        
        result=rbind(PhenotypesA, PhenotypesB)
        write.csv(result,paste0(DistOutputPath, TreeString, "_SiderophoreComparisons",toString(i), ".csv"))
        ComparisonsList[[i]] = CheckDiffsGenomes
        i=i+1
        
        genomes_to_compare = genomes_to_compare[!(genomes_to_compare %in% result$IsolateID) ]
        
      }else{
        genomes_to_compare = genomes_to_compare[genomes_to_compare!=currentgenome]
      }
    }}
  
  
  
  # Staphylokinase comparisons
  #####################
  print("Staphylokinase")
  ComparisonsList = c()
  i=1
  
  if(nrow(NeedToCompareStaphylokinase) > 0 ){
    genomes_to_compare = unique((NeedToCompareStaphylokinase %>% arrange(-StaphylokinaseDistance))$Genome1)
    

    while(length(genomes_to_compare)>0){
      currentgenome = genomes_to_compare[1]
      DistSubset = AllDistances %>% filter(Genome1==currentgenome)
      
      A_Phenotype_Genomes = (DistSubset %>% filter( (StaphylokinaseDistance <= .2) & (CopheneticDist < 1e-5)))$Genome2
      B_Phenotype_Genomes = (DistSubset %>% filter( (StaphylokinaseDistance > .2) & (CopheneticDist < 1e-5)))$Genome2
      
  
      # check that genomes in two different groups at least have .1 distance 
      CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
      CheckDiffsGenomes = paste(CheckDiffs$Genome1, CheckDiffs$Genome2, sep="_")
      
      ChooseToExclude =  (CheckDiffs[which(CheckDiffs$StaphylokinaseDistance <= .2),])
      while(nrow(ChooseToExclude)>0){
        
        
        # if there's a kind of 'edge case' where there's some match of A vs. B genome without much phenotype distance (not >.2)
        # Remove the one with the maximum cophenetic distance from other genomes that's the greatest
        ListPotentialExcludeA = sapply(ChooseToExclude$Genome1, as.character)
        ListPotentialExcludeB = sapply(ChooseToExclude$Genome2, as.character)
        
        
        if((length(ListPotentialExcludeA) + length(ListPotentialExcludeB)) > 0){
          PotentialExcludes = c(ListPotentialExcludeA, ListPotentialExcludeB)
          if(length(B_Phenotype_Genomes) ==1){
            PotentialExcludes = ListPotentialExcludeA
          }
          if(length(A_Phenotype_Genomes)==1){
            PotentialExcludes = ListPotentialExcludeB
          }
          maxdist = 0
          for(p in PotentialExcludes){
            current = max( ( AllDistances %>% filter(Genome1==p))$CopheneticDist  ) 
            if(current > maxdist){
              remove=p
              maxdist=current
            }
            
          }
        }
        B_Phenotype_Genomes = B_Phenotype_Genomes[B_Phenotype_Genomes!=remove]
        A_Phenotype_Genomes = A_Phenotype_Genomes[A_Phenotype_Genomes!=remove]
        CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
        CheckDiffsGenomes = paste(CheckDiffs$Genome1, CheckDiffs$Genome2, sep="_")
        
        ChooseToExclude =  (CheckDiffs[which(CheckDiffs$StaphylokinaseDistance < .2),])
        
      }
      if(length(ComparisonsList)>0){
        LengthsIntersect = sapply(1:length(ComparisonsList),function(x) length(intersect(ComparisonsList[[x]], CheckDiffsGenomes)))
      }else{
        LengthsIntersect=c()
      }
      A_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% A_Phenotype_Genomes))
      B_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% B_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes))
      
      ACheck= ( (max(A_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(A_Phenotype_GenomeComparisons$StaphylokinaseDistance) <=.2 ) )
      BCheck= ( (max(B_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(B_Phenotype_GenomeComparisons$StaphylokinaseDistance) <=.2 ) )
      
      if(ACheck==T & BCheck==T & length(B_Phenotype_Genomes)>0 & length(A_Phenotype_Genomes)>0 & !any(LengthsIntersect>0)){
        PhenotypesA = TreePlot$data %>% select(IsolateID, staphylokinase) %>% filter(IsolateID %in% A_Phenotype_Genomes) 
        PhenotypesA$Phenotype="A"
        
        PhenotypesB = TreePlot$data %>% select(IsolateID, staphylokinase) %>% filter(IsolateID %in% B_Phenotype_Genomes) 
        PhenotypesB$Phenotype="B"
        
        result=rbind(PhenotypesA, PhenotypesB)
        write.csv(result,paste0(DistOutputPath, TreeString, "_StaphylokinaseComparisons",toString(i), ".csv"))
        ComparisonsList[[i]] = CheckDiffsGenomes
        i=i+1
        
        genomes_to_compare = genomes_to_compare[!(genomes_to_compare %in% result$IsolateID) ]
      } else{
        genomes_to_compare = genomes_to_compare[genomes_to_compare!=currentgenome]
      }
      
    }}
  
  
  # Staphyloxanthin comparisons
  ###########################
  print("Staphyloxanthin")
  
  ComparisonsList = c()
  i=1
  
  if(nrow(NeedToCompareStaphyloxanthin) > 0 ){
    genomes_to_compare = unique((NeedToCompareStaphyloxanthin %>% arrange(-StaphyloxanthinDistance))$Genome1)

    while(length(genomes_to_compare)>0){
      currentgenome = genomes_to_compare[1]
      DistSubset = AllDistances %>% filter(Genome1==currentgenome)
      
      A_Phenotype_Genomes = (DistSubset %>% filter( (StaphyloxanthinDistance <= .2) & (CopheneticDist < 1e-5)))$Genome2
      B_Phenotype_Genomes = (DistSubset %>% filter( (StaphyloxanthinDistance > .2) & (CopheneticDist < 1e-5)))$Genome2
      

      # check that genomes in two different groups at least have .1 distance 
      CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
      CheckDiffsGenomes = paste(CheckDiffs$Genome1, CheckDiffs$Genome2, sep="_")
      ChooseToExclude =  (CheckDiffs[which(CheckDiffs$StaphyloxanthinDistance <= .2),])
      while(nrow(ChooseToExclude)>0){
        
        
        # if there's a kind of 'edge case' where there's some match of A vs. B genome without much phenotype distance (not >.2)
        # Remove the one with the maximum cophenetic distance from other genomes that's the greatest
        ListPotentialExcludeA = sapply(ChooseToExclude$Genome1, as.character)
        ListPotentialExcludeB = sapply(ChooseToExclude$Genome2, as.character)
        
        
        if((length(ListPotentialExcludeA) + length(ListPotentialExcludeB)) > 0 ){
          PotentialExcludes = c(ListPotentialExcludeA, ListPotentialExcludeB)
          if(length(B_Phenotype_Genomes) ==1){
            PotentialExcludes = ListPotentialExcludeA
          }
          if(length(A_Phenotype_Genomes)==1){
            PotentialExcludes = ListPotentialExcludeB
          }
          maxdist = 0
          for(p in PotentialExcludes){
            current = max( ( AllDistances %>% filter(Genome1==p))$CopheneticDist  ) 
            if(current > maxdist){
              remove=p
              maxdist=current
            }
            
          }
        }
        B_Phenotype_Genomes = B_Phenotype_Genomes[B_Phenotype_Genomes!=remove]
        A_Phenotype_Genomes = A_Phenotype_Genomes[A_Phenotype_Genomes!=remove]
        CheckDiffs = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes)) 
        CheckDiffsGenomes = paste(CheckDiffs$Genome1,CheckDiffs$Genome2, sep="_")
        ChooseToExclude =  (CheckDiffs[which(CheckDiffs$StaphyloxanthinDistance < .2),])
        
      }
      A_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% A_Phenotype_Genomes) & (Genome2 %in% A_Phenotype_Genomes))
      B_Phenotype_GenomeComparisons = AllDistances %>% filter((Genome1 %in% B_Phenotype_Genomes) & (Genome2 %in% B_Phenotype_Genomes))
      
      ACheck= ( (max(A_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(A_Phenotype_GenomeComparisons$StaphyloxanthinDistance) <=2 ) )
      BCheck= ( (max(B_Phenotype_GenomeComparisons$CopheneticDist) <1e-5) & (max(B_Phenotype_GenomeComparisons$StaphyloxanthinDistance) <=2 ) )
      
      
      if(length(ComparisonsList)>0){
        LengthsIntersect = sapply(1:length(ComparisonsList),function(x) length(intersect(ComparisonsList[[x]], CheckDiffsGenomes)))
      }else{
        LengthsIntersect=c()
      }
      
      if(ACheck==T & BCheck==T & length(B_Phenotype_Genomes)>0 & length(A_Phenotype_Genomes)>0 & !any(LengthsIntersect>0)){
        PhenotypesA = TreePlot$data %>% select(IsolateID, staphyloxanthin) %>% filter(IsolateID %in% A_Phenotype_Genomes) 
        PhenotypesA$Phenotype="A"
        
        PhenotypesB = TreePlot$data %>% select(IsolateID, staphyloxanthin) %>% filter(IsolateID %in% B_Phenotype_Genomes) 
        PhenotypesB$Phenotype="B"
        
        result=rbind(PhenotypesA, PhenotypesB)
        write.csv(result,paste0(DistOutputPath, TreeString, "_StaphyloxanthinComparisons",toString(i), ".csv"))
        ComparisonsList[[i]] = CheckDiffsGenomes
        i=i+1
        
        genomes_to_compare = genomes_to_compare[!(genomes_to_compare %in% result$IsolateID) ]
        
      }else{
        genomes_to_compare = genomes_to_compare[genomes_to_compare!=currentgenome]
      }
    }}
  
  
  
    # Visually, for these two it's obvious that we should compare SA1085/SA1081 vs. SA1086/SA1082
    # looking at NeedToCompareBiofilm, if we start with the first in order (1081):
    # Look at the AllDistances matrix comparison of 1081 to others:
    
    # genomes_to_compare = unique(NeedToCompareBiofilm$Genome1) = SA1081, SA1082, SA1085, SA1086
    
    # iterate through genomes_to_compare. 
    # While length(genomes_to_compare) > 0 :
    # currentgenome=genomes_to_compare[1]
    # (starting with SA1081 in this case)
    # Subset AllDistances %>% filter(Genome1 == "SA1081")
    # Identify genomes with <1e-5 & BiofilmDistance > .2
    
    # Set A_Phenotype_Genomes = genomes with <1e-5 & BiofilmDistance <= .2 as compared to SA1081 
    # Set B_Phenotype_Genomes = genomes with <1e-5 & BiofilmDistance > .2 as compared to SA1081 (including SA1081)
    
    # Make one A_Phenotypes subset of comparisons such that BOTH Genome1 and Genome2 are in A_Phenotype_Genomes
    # Make one B_Phenotypes subset of comparisons such that BOTH Genome1 and Genome2 are in B_Phenotypes_Genomes
    
    # set i=1
    # Set ACheck = whether no Biofilm distance in A_Phenotypes >.2 and no Cophenetic distance in A_Phenotypes >=1e-5
    # Set BCheck = whether no Biofilm distance in B_Phenotypes >.2 and no Cophenetic distance in B_Phenotypes >=1e-5
    # If ACheck==T & BCheck == T :
    #       Make a dataframe with rows c(IsolateID, biofilm, "PhenotypeB") for all the genomes in B 
    #       Make a dataframe with rows c(IsolateID, biofilm, "PhenotypeA") for all the genomes in A 
    #       Bind those dataframes together
    #       Then, save them as like TreeString_Biofilm_Comparison<i>.csv
    #       genomes_to_compare =genomes_to_compare[!(genomes_to_compare %in% c(A_Phenotypes$IsolateID, B_Phenotypes$IsolateID))]
    #       i=i+1
  
    
  }
  

  
  ##################################################################
  # Identify if any of the pairings have <1-e5 phylogenetic distance
  ##################################################################
  
  # colorpalettedf = colormap %>% filter(week %in%  TreePlot$data$week_collected)
  # TreePlot = TreePlot+ geom_tippoint(shape=15,size=5, aes(color=factor(week_collected))) +theme(legend.position="none") + scale_color_manual(values=colorpalettedf$colors)
  # genomesIncluded = TreeInput$tip.label
  # genomesIncluded = genomesIncluded[genomesIncluded!=RefGenome]
  # heightsave=6*(5/length(genomesIncluded))
  # 
  # TreePlot$data$staphylokinase[is.na(TreePlot$data$staphylokinase)] <- as.numeric(1.0)
  # TreePlot$data$staphyloxanthin[is.na(TreePlot$data$staphyloxanthin)] <- as.numeric(1.0)
  # TreePlot$data$biofilm[is.na(TreePlot$data$biofilm)] <- as.numeric(1.0)
  # TreePlot$data$siderophore[is.na(TreePlot$data$siderophore)] <- as.numeric(1.0)
  # 
  # TreePlot = TreePlot + geom_treescale(width=1e-6) + new_scale_fill()+ geom_fruit(aes(x=(staphyloxanthin)),fill="#B8860B",axis.params=list(add.axis=TRUE, axis="x"), geom=geom_bar, stat="identity", orientation="y", offset=.2,pwidth = .2)+
  #   new_scale_fill()+ geom_fruit(aes(x=(staphylokinase)),fill="#2D9D92", geom=geom_bar, stat="identity", orientation="y",axis.params=list(add.axis=TRUE, axis="x"),offset=0,pwidth = .2)+
  # new_scale_fill()+ geom_fruit(aes(x=(biofilm)),fill="#83B44B", geom=geom_bar, stat="identity", orientation="y",axis.params=list(add.axis=TRUE, axis="x"),offset=0,pwidth = .2)+
  # new_scale_fill()+ geom_fruit(aes(x=(siderophore)),fill="#8B0000", geom=geom_bar, stat="identity", orientation="y",axis.params=list(add.axis=TRUE, axis="x"),offset=0,pwidth = .2)
  # ggsave(TreePlot, file=paste0("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/ClonalFrameTreesAnnotated/", TreeString, "_tree.pdf"),height=heightsave, width=14.1)
  # 

  
  ###########################################
  ## Phage and plasmid presence by tree tip:
  ###########################################

  # plasmidContentsPatient = plasmidContents %>% filter(IsolateID %in% genomesIncluded )
  # plasmidContentsPatient = plasmidContentsPatient %>% select_if(function(x) any(x!=0))
  # idorder=plasmidContentsPatient$IsolateID
  # plasmidContentsPatient = plasmidContentsPatient %>% select_if(function(x) any(x==0)) # differentially present/absent plasmids only
  # plasmidContentsPatient$IsolateID = idorder
  # plasmidContentsPatient$type = "Plasmid"
  # plasmidcontentsmelt= plasmidContentsPatient %>% melt(id.vars=c("IsolateID","type"))
  #   
  # phageContentsPatient = phageContents %>% filter(IsolateID %in% genomesIncluded )
  # phageContentsPatient = phageContentsPatient %>% select_if(function(x) any(x!=0))
  # isolateorderphage = phageContentsPatient$IsolateID
  # phageContentsPatient = phageContentsPatient %>% select_if(function(x) any(x==0)) 
  # phageContentsPatient$IsolateID = isolateorderphage
  # phageContentsPatient$type="Phage"
  # phageContentsMelt = phageContentsPatient %>% melt(id.vars=c("IsolateID","type"))
  # 
  # 
  # if( (ncol(phageContentsMelt)) > 2 & (ncol(plasmidcontentsmelt) > 2 )){
  #   phageContentsMelt$value[phageContentsMelt$value==1] <- 2
  #   HGT_Info = rbind(plasmidcontentsmelt, phageContentsMelt)
  #   starshapeobj=c(11,15)
  #   fillvals=c("white","darkblue","#024B30")
  #   colvals = c("white","white","#024B30")
  #   
  # }
  # if((ncol(phageContentsMelt) > 2) & (ncol(plasmidcontentsmelt)==2) ){
  #   phageContentsMelt$value[phageContentsMelt$value==1] <- 2
  #   HGT_Info=phageContentsMelt
  #   starshapeobj=c(11)
  #   fillvals=c("white","#024B30")
  #   colvals = c("white","#024B30")
  # }
  # if((ncol(phageContentsMelt)==2) & (ncol(plasmidcontentsmelt)>2 )){
  #   HGT_Info=plasmidcontentsmelt
  #   fillvals=c("white","darkblue")
  #   colvals = c("white", "white")
  #   starshapeobj=c(15)
  #   
  # }
  # if( (ncol(phageContentsMelt)==2) & (ncol(plasmidcontentsmelt)==2)){
  #   HGT_Info=data.frame()
  #   print(paste0("No variably present/absent phages or plasmids in ",TreeString ))
  # }
  # 
  # if(nrow(HGT_Info) > 0){
  #   
  #   yorder = (TreePlot$data %>% arrange(y) %>% filter(grepl(label, pattern="SA") & (label!=RefGenome)))$label
  #   
  #   
  #   widthobj=10*(length(unique(HGT_Info$variable))/ 9)
  #   HGTplot = ggplot(HGT_Info, aes(y=IsolateID, x=variable, color=factor(value), fill=factor(value), starshape=factor(type)))+ geom_star(size=1,size=1) +scale_color_manual(values=colvals)+ scale_starshape_manual(values=starshapeobj)+
  #    scale_fill_manual(values=fillvals) + theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), axis.text.y=element_text(size=2),
  #                                               legend.position="none", axis.line.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
  #                                               panel.background = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())+scale_x_discrete(position = "top") #+ coord_fixed(ratio=1)
  #   
  #   ThingOrder = sapply(unique(HGTplot$data$variable),as.character)
  #   HGTplot$data$IsolateID = factor(HGTplot$data$IsolateID, levels=c(yorder, setdiff(megadf$IsolateID, yorder)))
  #   HGTplot$data$variable = factor(HGTplot$data$variable, levels=c(ThingOrder, setdiff(allplasmid_phages, ThingOrder)))
  #   
  #   
  #   HGTplot = HGTplot + scale_y_discrete(limits= c(yorder, setdiff(megadf$IsolateID, yorder)), drop=FALSE) + scale_x_discrete(limits= c(ThingOrder, setdiff(allplasmid_phages, ThingOrder)), drop=FALSE)
  #   
  #   
  #   ggsave(HGTplot, file=paste0("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/ClonalFrameTreesAnnotated/", TreeString, "_HGTplot.pdf"),height=12.5,width=7)

  #}





