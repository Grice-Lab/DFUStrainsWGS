# Amy Campbell
# 11/2020
# Zooming in on the 'clades' as defined by <500 snp distance which map to cc8, cc1, cc5 (most common)

# cc5 look below 357
# cc8 look below 307
# cc1 look below 260

SubjTime = function(row){
  if(is.na(row["DORN"])){
    return(NA)
  }else if(is.na(row["patient_id"])){
    return(row["DORN"])
  }else{
    return(paste(row["patient_id"], stringr::str_replace(row["visit"], " ", ""), sep="_"))
  }}


setwd('/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/')
library("treeio")
library("ggtree")
library("dplyr")
library("ggplot2")
library("dplyr")
library("ape")
library("tibble")

Isolates = "data/DFU_Staph_aureus_isolates.csv"
ARM_Phenotypes = read.csv("data/Phenotypes_01.15.20.csv")
StaphIsolateDORNs=read.csv("data/DFU_Staph_aureus_isolates.csv")
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




apeplot <- ggtree::ggtree(apeRooted, size=.25, layout = "rectangular") + geom_tiplab(size=2)
apeplotbackup = apeplot
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


fullPlotData = Plot_Data[c("label", "DORN", "patient_id", "visit", "TipLABEL", "SubjTime", "HealedYorN", "WeekHealed", "HealingTime")]

# CC5
#####
# cc5 look below 357
ggtree::getSubtree(apeRooted, 357)

CC5_ape_subset = ape::extract.clade(apeRooted, 357)
CC5_Tree = ggtree(CC5_ape_subset, size=.25, layout = "rectangular") + geom_tiplab(size=2)

CC5_Tree$data = CC5_Tree$data %>% left_join(fullPlotData,by="label") 
CC5_Tree$data$label1 = CC5_Tree$data$label 
CC5_Tree$data$label = CC5_Tree$data$TipLABEL
CC5_Tree$data$HealedBy12  = if_else(CC5_Tree$data$WeekHealed <= 12, "Yes", "No")
ggtree(CC5_Tree$data, size=.25, layout = "rectangular") + geom_tippoint(aes(color=HealedBy12)) + geom_tiplab(size=3) + scale_color_manual(values=c(c("darkred", "#9FE2BF")))
length(unique(CC5_Tree$data$patient_id))

new <- CC5_Tree$data$DORN[!is.na(CC5_Tree$data$DORN)]

# CC8
#####
# cc8 look below 307
ggtree::getSubtree(apeRooted, 307)

CC8_ape_subset = ape::extract.clade(apeRooted, 307)
CC8_Tree = ggtree(CC8_ape_subset, size=.25, layout = "rectangular") + geom_tiplab(size=2)

CC8_Tree$data = CC8_Tree$data %>% left_join(fullPlotData,by="label") 
CC8_Tree$data$label1 = CC8_Tree$data$label 
CC8_Tree$data$label = CC8_Tree$data$TipLABEL
CC8_Tree$data$HealedBy12  = if_else(CC8_Tree$data$WeekHealed <= 12, "Yes", "No")
labeled_healing = ggtree(CC8_Tree$data, size=.25, layout = "rectangular") + geom_tippoint(aes(color=HealedBy12)) + geom_tiplab(size=3) + scale_color_manual(values=c(c("darkred", "#9FE2BF")))

length(unique(CC8_Tree$data$patient_id))

na.omit(CC8_Tree$data$label)

ggsave(labeled_healing, file="CC8Tree_Labeled_healing.pdf", width=10, height=8)
labeled_healing$data$label= CC8_Tree$data$DORN
ggsave(labeled_healing, file="CC8Tree_Labeled_DORN.pdf", width=10, height=8)
CC8_DORNS = CC8_Tree$data$DORN[!is.na(CC8_Tree$data$DORN)]
CC8_DORNS
write.csv(CC8_DORNS, "CC8_Genomes.txt", row.names=F, col.names=F)

# cc1
#####
# cc1 look below 307
ggtree::getSubtree(apeRooted, 260)

cc1_ape_subset = ape::extract.clade(apeRooted, 260)
cc1_Tree = ggtree(cc1_ape_subset, size=.25, layout = "rectangular") + geom_tiplab(size=2)

cc1_Tree$data = cc1_Tree$data %>% left_join(fullPlotData,by="label") 
cc1_Tree$data$label1 = cc1_Tree$data$label 
cc1_Tree$data$label = cc1_Tree$data$TipLABEL
cc1_Tree$data$HealedBy12  = if_else(cc1_Tree$data$WeekHealed <= 12, "Yes", "No")
ggtree(cc1_Tree$data, size=.25, layout = "rectangular") + geom_tippoint(aes(color=HealedBy12)) + geom_tiplab(size=3) + scale_color_manual(values=c(c("darkred", "#9FE2BF")))
length(unique(cc1_Tree$data$patient_id))




commonsizedplasmid = CC5_Tree$data %>% filter(DORN %in% c("DORN2144",
                                                          "DORN2187",
                                                          "DORN2083",
                                                          "DORN2089",
                                                          "DORN2123",
                                                          "DORN1289",
                                                          "DORN1285"))

UP_1405 = c("DORN1881",
            "DORN1702",
            "DORN1740",
            "DORN1761",
            "DORN1779",
            "DORN315",
            "DORN2083",
            "DORN2089",
            "DORN2123",
            "DORN2144",
            "DORN2187",
            "DORN807")
CC5_Tree$data %>% filter((DORN) %in% UP_1405)

CC5_Tree$data %>% mutate(UP1405_plasmid = if_else((DORN) %in% UP_1405, TRUE,FALSE ))
cc5tree = ggtree(CC5_Tree$data, size=.25, layout = "rectangular") + geom_tippoint(aes(color=HealedBy12)) + geom_tiplab(size=3) + scale_color_manual(values=c(c("darkred", "#9FE2BF"))) 
ggsave(cc5tree, file="cc5_TreeHealing.pdf", height=20, width=15)
cc5tree$data$label = cc5tree$data$DORN
ggsave(cc5tree, file="cc5_TreeHealingDORN.pdf", height=20, width=15)
