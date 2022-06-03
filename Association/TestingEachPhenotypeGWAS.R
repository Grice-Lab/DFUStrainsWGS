# Amy Campbell
# 06-01-2022
# Retesting GWAS systematically comparing each with no subsetting vs. subsetting by CC, each phenotype
library(dplyr)
library(reshape2)
library(ggplot2)
library(KernSmooth)

AllPhenotypes=read.csv("data/staphyloxanthin_paper_data.csv")
AllPhenotypes$DORN = paste0("DORN", AllPhenotypes$DORN)
CCmapping= read.csv("data/Phylogeny2022Data/CCMapPlotting.csv")
KinasePresence = read.csv("data/staphyloxanthin_paper_UpdatedSakClassifications.csv")
KinasePresent = (KinasePresence %>% filter(SakUpdated == "yes"))$DORN
# Staphyloxanthin production
############################

# not using DORN429 or DORN1176
AllPhenotypes = AllPhenotypes %>% filter(!(DORN %in% c("DORN429", "DORN1176")))
dim(AllPhenotypes)
setdiff(CCmapping$DORN, AllPhenotypes$DORN)


AllPhenotypes = AllPhenotypes %>% left_join(CCmapping, by="DORN")


AllPhenotypesZeroOne = AllPhenotypes
AllPhenotypesZeroOne = AllPhenotypesZeroOne %>% select(staphyloxanthin,staphylokinase, siderophore, biofilm, DORN, patient, CCLabel)


# Zero-one normalize all four phenotypes 
#########################################
AllPhenotypesZeroOne$staphyloxanthin = ((AllPhenotypesZeroOne$staphyloxanthin  - min(AllPhenotypesZeroOne$staphyloxanthin))/(max(AllPhenotypesZeroOne$staphyloxanthin) - min(AllPhenotypesZeroOne$staphyloxanthin)))
AllPhenotypesZeroOne$staphylokinase = ((AllPhenotypesZeroOne$staphylokinase  - min(AllPhenotypesZeroOne$staphylokinase, na.rm=T))/(max(AllPhenotypesZeroOne$staphylokinase, na.rm = T) - min(AllPhenotypesZeroOne$staphylokinase,na.rm = T)))
AllPhenotypesZeroOne$biofilm = ((AllPhenotypesZeroOne$biofilm  - min(AllPhenotypesZeroOne$biofilm))/(max(AllPhenotypesZeroOne$biofilm) - min(AllPhenotypesZeroOne$biofilm)))
AllPhenotypesZeroOne$siderophore = ((AllPhenotypesZeroOne$siderophore  - min(AllPhenotypesZeroOne$siderophore))/(max(AllPhenotypesZeroOne$siderophore) - min(AllPhenotypesZeroOne$siderophore)))


NormedXanthin = (AllPhenotypesZeroOne %>% select(DORN, staphyloxanthin))
NormedXanthin$Path = paste0("/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/", NormedXanthin$DORN)
NormedXanthin$Path = paste0(NormedXanthin$Path, "_Final.fasta")
colnames(NormedXanthin) = c("ID", "Phenotype", "Path")
write.table(NormedXanthin, file="mappings/TestZeroOneStaphyloxanthin.txt", sep="\t", row.names=F, quote=F)

NormedSiderophore= (AllPhenotypesZeroOne %>% select(DORN, siderophore))
NormedSiderophore$Path = paste0("/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/", NormedSiderophore$DORN)
NormedSiderophore$Path = paste0(NormedSiderophore$Path, "_Final.fasta")
write.table(NormedSiderophore, file="mappings/TestZeroOneSiderophore.txt", sep="\t", row.names=F, quote=F)

NormedStaphylokinase = AllPhenotypesZeroOne %>% filter(DORN %in% KinasePresent) %>% select(DORN, staphylokinase)
NormedStaphylokinase$Path = paste0("/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/", NormedStaphylokinase$DORN)
NormedStaphylokinase$Path = paste0(NormedStaphylokinase$Path, "_Final.fasta")
NormedStaphylokinase= NormedStaphylokinase %>% filter(!is.na(staphylokinase))
write.table(NormedStaphylokinase, file="mappings/TestZeroOneStaphylokinase.txt", sep="\t", row.names=F, quote=F)


NormedBiofilm= (AllPhenotypesZeroOne %>% select(DORN, biofilm))
NormedBiofilm$Path = paste0("/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/", NormedSiderophore$DORN)
NormedBiofilm$Path = paste0(NormedBiofilm$Path, "_Final.fasta")
write.table(NormedBiofilm, file="mappings/TestZeroOneBiofilm.txt", sep="\t", row.names=F, quote=F)


AllPhenotypesUnique = AllPhenotypes %>% select(patient, CCLabel) %>% unique()
plots=rep(0, 70)
i=1
for(j in 1:nrow(AllPhenotypesUnique)){
  
  Patient = (AllPhenotypesUnique[j, "patient"])
  CC = (AllPhenotypesUnique[j, "CCLabel"])
  Subset = AllPhenotypesZeroOne %>% filter(CCLabel==CC & patient==Patient)
  if(nrow(Subset)>1){
    SubsetMelt = Subset %>% select(patient, CCLabel, DORN, staphylokinase, staphyloxanthin, siderophore, biofilm) %>% reshape2::melt(id.vars=c("patient","CCLabel","DORN"))
    ggsave(ggplot(SubsetMelt, aes(x=DORN, y=value)) + geom_point() + ylim(0,1)+ facet_grid(.~ variable) + theme(axis.text.x=element_text(angle=270)), file=paste0("data/Phylogeny2022Data/", Patient, CC, ".pdf"))
    i=i+1
  }
}


# Test
########################################


biofilmtest = sort(AllPhenotypesZeroOne$biofilm)
siderophoretest = sort(AllPhenotypesZeroOne$siderophore)
saktest = sort(NormedStaphylokinase$staphylokinase)
staphyloxanthintest = sort(AllPhenotypesZeroOne$staphyloxanthin)

ks_biofilm = KernSmooth::bkde(biofilmtest, kernel = "normal", canonical = FALSE, .025,
     gridsize = 401,c(0,1), truncate = TRUE)
ks_siderophore = KernSmooth::bkde(siderophoretest, kernel = "normal", canonical = FALSE, .025,
                              gridsize = 401,c(0,1), truncate = TRUE)
ks_staphylokinase = KernSmooth::bkde(saktest, kernel = "normal", canonical = FALSE, .025,
                                  gridsize = 401,c(0,1), truncate = TRUE)

ks_staphyloxanthin = KernSmooth::bkde(staphyloxanthintest, kernel = "normal", canonical = FALSE, .025,
                                      gridsize = 401,c(0,1), truncate = TRUE)

biofilmkde = ggplot() + geom_bar(data=AllPhenotypesZeroOne, aes(x=biofilm),stat="count") + geom_line(data=data.frame(ks_biofilm), aes(x=x, y=y), color="red") + ggtitle("KDE Biofilm")
sakKDE= ggplot() +geom_bar(data=NormedStaphylokinase, aes(x=staphylokinase),stat="count")+ geom_line(data=data.frame(ks_staphylokinase), aes(x=x, y=y), color="red") + ggtitle("KDE Staphylokinase") # geom_bar(aes(y=Sales,fill="Sales"), width=.7, stat="identity")
sidKDE = ggplot() +geom_bar(data=AllPhenotypesZeroOne, aes(x=siderophore),stat="count")+ geom_line(data=data.frame(ks_siderophore), aes(x=x, y=y),color="red") + ggtitle("KDE Siderophore") # geom_bar(aes(y=Sales,fill="Sales"), width=.7, stat="identity")
xanthinKDE = ggplot() +geom_bar(data=AllPhenotypesZeroOne, aes(x=staphyloxanthin),stat="count")+ geom_line(data=data.frame(ks_staphyloxanthin), aes(x=x, y=y),color="red") + ggtitle("KDE Staphyloxanthin") # geom_bar(aes(y=Sales,fill="Sales"), width=.7, stat="identity")

gridExtra::grid.arrange(biofilmkde,sakKDE,  sidKDE, xanthinKDE)

IDminima = function(vector_x, vector_y){
  minima = c()
  for(j in 1:length(vector_y)){
    if(j==1){
      minima = append(minima, vector_x[j])
    }else if(j==length(vector_y)){
      minima=append(minima, vector_x[j])
    }else{
      if( (vector_y[j-1] > vector_y[j]) & vector_y[j+1] > vector_y[j]){
        minima= append(minima, vector_x[j])
      }
      
    }
  }
  return(minima)
}

BiofilmCutoffs = IDminima(ks_biofilm$x, ks_biofilm$y)
SiderophoreCutoffs = IDminima(ks_siderophore$x, ks_siderophore$y)
XanthinCutoffs = IDminima(ks_staphyloxanthin$x, ks_staphyloxanthin$y)
KinaseCutoffs = IDminima(ks_staphylokinase$x, ks_staphylokinase$y)


AllPhenotypesZeroOne = AllPhenotypesZeroOne %>% mutate(biofilm_cluster = case_when( (biofilm >= BiofilmCutoffs[1]) & (biofilm < BiofilmCutoffs[2]) ~ "1",
                                                                                    (biofilm >= BiofilmCutoffs[2]) & (biofilm < BiofilmCutoffs[3]) ~ "2",
                                                                                    (biofilm >= BiofilmCutoffs[3]) & (biofilm < BiofilmCutoffs[4]) ~ "3",
                                                                                    (biofilm >= BiofilmCutoffs[4]) & (biofilm < BiofilmCutoffs[5]) ~ "4",
                                                                                    (biofilm >= BiofilmCutoffs[5]) & (biofilm < BiofilmCutoffs[6]) ~ "5", 
                                                                                    (biofilm >= BiofilmCutoffs[6]) & (biofilm < BiofilmCutoffs[7]) ~ "6", 
                                                                                    (biofilm >= BiofilmCutoffs[7]) & (biofilm < BiofilmCutoffs[8]) ~ "7",
                                                                                    (biofilm >= BiofilmCutoffs[8]) & (biofilm <= BiofilmCutoffs[9]) ~ "8"))

AllPhenotypesZeroOne = AllPhenotypesZeroOne %>% mutate(siderophore_cluster = case_when( (siderophore >= SiderophoreCutoffs[1]) & (siderophore < SiderophoreCutoffs[2]) ~ "1",
                                                                                    (siderophore >= SiderophoreCutoffs[2]) & (siderophore < SiderophoreCutoffs[3]) ~ "2",
                                                                                    (siderophore >= SiderophoreCutoffs[3]) & (siderophore < SiderophoreCutoffs[4]) ~ "3",
                                                                                    (siderophore >= SiderophoreCutoffs[4]) & (siderophore < SiderophoreCutoffs[5]) ~ "4",
                                                                                    (siderophore >= SiderophoreCutoffs[5]) & (siderophore < SiderophoreCutoffs[6]) ~ "5", 
                                                                                    (siderophore >= SiderophoreCutoffs[6]) & (siderophore < SiderophoreCutoffs[7]) ~ "6", 
                                                                                    (siderophore >= SiderophoreCutoffs[7]) & (siderophore < SiderophoreCutoffs[8]) ~ "7",
                                                                                    (siderophore >= SiderophoreCutoffs[8]) & (siderophore <= SiderophoreCutoffs[9]) ~ "8"))


AllPhenotypesZeroOne = AllPhenotypesZeroOne %>% mutate(kinase_cluster = case_when( (staphylokinase >= KinaseCutoffs[1]) & (staphylokinase < KinaseCutoffs[2]) ~ "1",
                                                                                        (staphylokinase >= KinaseCutoffs[2]) & (staphylokinase < KinaseCutoffs[3]) ~ "2",
                                                                                        (staphylokinase >= KinaseCutoffs[3]) & (staphylokinase < KinaseCutoffs[4]) ~ "3",
                                                                                        (staphylokinase >= KinaseCutoffs[4]) & (staphylokinase < KinaseCutoffs[5]) ~ "4",
                                                                                        (staphylokinase >= KinaseCutoffs[5]) & (staphylokinase < KinaseCutoffs[6]) ~ "5", 
                                                                                        (staphylokinase >= KinaseCutoffs[6]) & (staphylokinase < KinaseCutoffs[7]) ~ "6", 
                                                                                        (staphylokinase >= KinaseCutoffs[7]) & (staphylokinase <= KinaseCutoffs[8]) ~ "7",
                                                                                         ))

AllPhenotypesZeroOne = AllPhenotypesZeroOne %>% mutate(xanthin_cluster = case_when( (staphyloxanthin >= XanthinCutoffs[1]) & (staphyloxanthin < XanthinCutoffs[2]) ~ "1",
                                                                                   (staphyloxanthin >= XanthinCutoffs[2]) & (staphyloxanthin < XanthinCutoffs[3]) ~ "2",
                                                                                   (staphyloxanthin >= XanthinCutoffs[3]) & (staphyloxanthin < XanthinCutoffs[4]) ~ "3",
                                                                                   (staphyloxanthin >= XanthinCutoffs[4]) & (staphyloxanthin < XanthinCutoffs[5]) ~ "4",
                                                                                   (staphyloxanthin >= XanthinCutoffs[5]) & (staphyloxanthin < XanthinCutoffs[6]) ~ "5", 
                                                                                   (staphyloxanthin >= XanthinCutoffs[6]) & (staphyloxanthin < XanthinCutoffs[7]) ~ "6", 
                                                                                   (staphyloxanthin >= XanthinCutoffs[7]) & (staphyloxanthin <= XanthinCutoffs[8]) ~ "7",
))

AllPhenotypesUnique = AllPhenotypesZeroOne %>% select(patient, CCLabel) %>% unique()
plots2=rep(0, 70)
i=1
for(j in 1:nrow(AllPhenotypesUnique)){
  
  Patient = (AllPhenotypesUnique[j, "patient"])
  CC = (AllPhenotypesUnique[j, "CCLabel"])
  Subset = AllPhenotypesZeroOne %>% filter(CCLabel==CC & patient==Patient)
  if(nrow(Subset)>1){
    SubsetMelt = Subset %>% select(patient, CCLabel, DORN, biofilm, biofilm_cluster) #%>% reshape2::melt(id.vars=c("patient","CCLabel","DORN", "biofilm_cluster",)) %>% filter()
    ggsave(ggplot(SubsetMelt, aes(x=DORN, y=biofilm, color=biofilm_cluster)) + geom_point() + ylim(0,1) + theme(axis.text.x=element_text(angle=270)), file=paste0("data/Phylogeny2022Data/MultipleGenomes/Biofilm", Patient, CC, ".pdf"))
    i=i+1
  }
}

dim(AllPhenotypesZeroOne %>% select(patient, CCLabel, biofilm_cluster) %>% unique())
dim(AllPhenotypesZeroOne %>% select(patient, CCLabel, siderophore_cluster) %>% unique())
dim(AllPhenotypesZeroOne %>% select(patient, CCLabel, kinase_cluster) %>% unique())
dim(AllPhenotypesZeroOne %>% select(patient, CCLabel, xanthin_cluster) %>% unique())


i=1
AllPhenotypesUniqueBiofilm = AllPhenotypesZeroOne %>% select(CCLabel, patient, biofilm_cluster) %>% unique()
UniqueIsolatesBiofilm =c()
RequireChoice = c()
RequireChoiceDORNs=c()
for(j in 1:nrow(AllPhenotypesUniqueBiofilm)){
  
  Patient = (AllPhenotypesUniqueBiofilm[j, "patient"])
  CC = (AllPhenotypesUniqueBiofilm[j, "CCLabel"])
  cluster = (AllPhenotypesUniqueBiofilm[j, "biofilm_cluster"])
  
  Subset = AllPhenotypesZeroOne %>% filter(CCLabel==CC & patient==Patient & biofilm_cluster==cluster)
  if(nrow(Subset)==1){
    UniqueIsolatesBiofilm = append(UniqueIsolatesBiofilm, Subset$DORN)
  }else{
    stringSubset = (paste(CC, Patient, cluster,sep="_"))
    list(Subset$DORN)
    SubsetDORNs = paste(Subset$DORN,collapse = "_")
    print(SubsetDORNs)
    RequireChoice= append(RequireChoice, stringSubset)
    RequireChoiceDORNs  = append(RequireChoiceDORNs, SubsetDORNs)
  }
  

}
BiofilmSelectDF = data.frame(Combo=RequireChoice, DORNs = RequireChoiceDORNs)

# Write a csv that i'll then manually add the chosen DORNs for each to  
write.csv(BiofilmSelectDF, file="data/Phylogeny2022Data/MultipleGenomes/BiofilmSelect.csv")
BiofilmSelectedDF =read.csv("data/Phylogeny2022Data/BiofilmSelected.csv")
BiofilmSelectedDF$Selected
UniqueIsolatesBiofilm
BiofilmSelectedList = append(BiofilmSelectedDF$Selected,UniqueIsolatesBiofilm )
NormedBiofilmSubset = NormedBiofilm %>% filter(DORN %in% BiofilmSelectedList)
write.table(NormedBiofilmSubset,file="mappings/BiofilmZeroOneSubsetKDE.txt", sep="\t", row.names=F, quote=F)
