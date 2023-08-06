library(dplyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)

randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#F6BE00",
                "#33A02C","#E6F5C9")
randpalette14 = randpalette18[-c(3,6,8,15)]

SplitKegg = function(x, column){
  # Column should be 1(ortholog name) or 2(Kegg ID)
  if(grepl(",", x)==T){
    col1 = str_split(x, ",")[[1]][1]
    col2 = str_split(x, ",")[[1]][2]
  } else{
    col1 = x
    col2 = ""
  }
  if(column==1){
    return(col1)
  }else{
    return(col2)
  }
}

Patient_Genome_Info = read.csv("~/Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")
Patient_Genome_Info$DORN = paste0("DORN", Patient_Genome_Info$Doern.lab.bank.)
patient_genome = Patient_Genome_Info %>% select(patient_id, DORN)

WeekHealed = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/staphyloxanthin_paper_data.csv") %>% select(patient, week_healed)
WeekHealed$HealedBy12 = if_else(WeekHealed$week_healed>12 | is.na(WeekHealed$week_healed), "No", "Yes")
HealedPatients = WeekHealed %>% filter(HealedBy12=="Yes")
UnhealedPatients = WeekHealed %>% filter(HealedBy12=="No")
Patient_Genome_Info$patient = Patient_Genome_Info$patient_id
Genomes_healing = Patient_Genome_Info %>% select(DORN, patient) %>% left_join(WeekHealed,  by="patient") %>% unique() %>% select(HealedBy12, DORN)


keggAnnotations = read.table("/Users/amycampbell/Documents/DataInputGithub/data/PanGenomeKegg.csv")

keggAnnotations$Gene = sapply(keggAnnotations$V1, function(x) SplitKegg(x, 1))
keggAnnotations$Kegg = sapply(keggAnnotations$V1, function(x) SplitKegg(x, 2))

# CCs to isolates
#################
CCMap = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Phylogeny2022Data/CCMapPlotting.csv")

################################################################
# Mapping annotations to gene presence/absence matrix from roary
################################################################ 
Blast2GoAnnotationFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoAnnotations_Original.txt", sep="\t",header=T,row.names=NULL) %>% select(GO.IDs, SeqName)
Blast2GoMapFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoMappingResults_Original.txt", sep="\t",header=T,row.names=NULL)# %>% select(GO.IDs, SeqName)
Blast2GoAnnotationFileNew = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoAnnotations_2149.txt", sep="\t",header=T,row.names=NULL) %>% select(GO.IDs, SeqName)

# SeqNames for the 'new' run of blast2go are the same as seqnames in most recent roary run
Blast2GoAnnotationFileNew$NewpangenomeID = Blast2GoAnnotationFileNew$SeqName
# But seqnames for the 'old' run of blast2go are the old roary run
Blast2GoAnnotationFile$OldPangenomeID = Blast2GoAnnotationFile$SeqName

# Map of pan-genome including SA2149 (new) to pan-genome excluding it(old)
OldToNewMap=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Map_New_To_Old_Pangenome.csv")

Blast2GoAnnotationFileNew = Blast2GoAnnotationFileNew %>% left_join(OldToNewMap, by="NewpangenomeID")
Blast2GoAnnotationFile = Blast2GoAnnotationFile %>% left_join(OldToNewMap, by="OldPangenomeID") %>% filter(!is.na(NewpangenomeID))

Blast2GoAnnotationFileNew = Blast2GoAnnotationFileNew %>% select(GO.IDs,NewpangenomeID)
Blast2GoAnnotationFile = Blast2GoAnnotationFile %>% select(GO.IDs, NewpangenomeID)
FullAnnotated  = rbind(Blast2GoAnnotationFileNew, Blast2GoAnnotationFile)

#Blast2GoMappingFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoAnnotations_Original.txt", sep="\t",header=T,row.names=NULL)
RoaryOutput = read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new.csv", sep=',')
RoaryOutput$NewpangenomeID = 1:nrow(RoaryOutput)

# mapping old to new in the roary gene_presence_absence file
# uncomment if you want to rerun this but it takes a little while
#################################################################
# For each item in OldToNewMap$NewpangenomeID, I want to figure out which row has an entry that == it and then stick it into a variable in that row
#for(id in OldToNewMap$NewpangenomeID){
  #listobj=apply(RoaryOutput, 2, function(x) which(x == id))
 # Index = (which(apply(RoaryOutput, 1, function(x) any(x==id))))
#  RoaryOutput[Index,"NewpangenomeID"] = id
#}
#RoaryOutput = RoaryOutput %>% left_join(OldToNewMap, by="NewpangenomeID")
#write.csv(RoaryOutput, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv")


RoaryComplete = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv")

RoaryOutputPresenceAbsence = RoaryComplete
RoaryOutputPresenceAbsence = RoaryOutputPresenceAbsence %>% select(Gene, OldPangenomeID, NewpangenomeID, colnames(RoaryOutputPresenceAbsence)[which(grepl("DORN", colnames(RoaryOutputPresenceAbsence)))])


RoaryOutputPresenceAbsence = RoaryOutputPresenceAbsence %>% mutate_at(vars(contains("DORN")), function(x) if_else(x=="", 0,1))

#####################################################
# Genes present on phages or plasmids we're tracking
####################################################
all_accessory = RoaryComplete
Binarized = all_accessory[15:(ncol(all_accessory)-2)] 
Binarized[Binarized==""] <-0
Binarized[Binarized!=0] <-1
all_accessory[15:(ncol(all_accessory)-2)]  <- Binarized %>% mutate_all(function(x) as.numeric(as.character(x)))

allcols = colnames(all_accessory)
colselect = c(allcols[grepl( "DORN",allcols)],"Gene")
all_accessory = all_accessory %>% select(colselect)

all_accessory$NumGenomes = rowSums(all_accessory[1:(ncol(all_accessory)-1)])

StrictCoreGenes = (all_accessory %>% filter(NumGenomes==220))$Gene
LessStrictCoreGenes = (all_accessory %>% filter(NumGenomes>=(220*.99)))$Gene
AccessoryGenes = (all_accessory %>% filter(NumGenomes<(220*.99)))$Gene

Genes_phage_plasmid = read.csv("~/Documents/DataInputGithub/data/IntraPatient/GenePresence_ByPhage_Plasmid_Updated.csv")
Genes_phage_plasmid$X.1=NULL
AccessoryGenePhagePlasmid = Genes_phage_plasmid%>% filter(X %in% AccessoryGenes)
AccessoryGenePhagePlasmid$NumOccurences = rowSums(AccessoryGenePhagePlasmid[,2:(ncol(AccessoryGenePhagePlasmid)-2)])
dim(AccessoryGenePhagePlasmid %>% filter(NumOccurences > 0))


CoregenePhagePlasmid =  Genes_phage_plasmid%>% filter(X %in% LessStrictCoreGenes)
CoregenePhagePlasmid$X.1=NULL
CoregenePhagePlasmid$NumOccurences = rowSums(CoregenePhagePlasmid[2:(ncol(CoregenePhagePlasmid)-2)])
dim(CoregenePhagePlasmid %>% filter(NumOccurences > 0))



######################################
# Ordination of accessory gene content 
######################################
RoaryOutputPresenceAbsenceForTransposition = RoaryOutputPresenceAbsence %>% select(colnames(RoaryOutputPresenceAbsence)[which(grepl("DORN", colnames(RoaryOutputPresenceAbsence)))])
RoaryOutputPresenceAbsenceForTranspositionTest = RoaryOutputPresenceAbsenceForTransposition

RoaryOutputPresenceAbsenceForTranspositionTest$Sums = rowSums(RoaryOutputPresenceAbsenceForTranspositionTest)
RoaryOutputPresenceAbsenceForTranspositionTest$Gene = RoaryOutputPresenceAbsence$Gene

AccessoryGenes = RoaryOutputPresenceAbsenceForTranspositionTest %>% filter(Sums <220) %>% filter(Sums>0)

AccessoryGenes$Sums = NULL
SaveGenes = AccessoryGenes$Gene
AccessoryGenes$Gene=NULL
GenomeNames = colnames(AccessoryGenes)

Presence_Absence_Accessory = data.frame(t(AccessoryGenes))
colnames(Presence_Absence_Accessory) = SaveGenes


JacDist = vegan::vegdist(Presence_Absence_Accessory,method="jaccard", binary=T )
PCOAresult = ape::pcoa(JacDist) 
PCOAresultVectors = data.frame(PCOAresult$vectors)
PCOAresultVectors$DORN = row.names(PCOAresultVectors)
PCOAresultVectors = PCOAresultVectors %>% left_join(CCMap, by="DORN")

Axis1_PctVar = round( (100* (PCOAresult$values$Relative_eig)[1]), 1)
Axis2_PctVar = round( (100* (PCOAresult$values$Relative_eig)[2]), 1)
Axis3_PctVar = round( (100* (PCOAresult$values$Relative_eig)[3]), 1)
Axis4_PctVar = round( (100* (PCOAresult$values$Relative_eig)[4]), 1)

PCOAresultVectors = PCOAresultVectors %>% left_join(Genomes_healing, by="DORN") 

PCOAresultVectors$Axis.1

PCoa_Plot_By_CC = ggplot(PCOAresultVectors, aes(x=Axis.1, y=Axis.2, color=factor(CCLabel), shape=factor(HealedBy12))) + geom_point(alpha=.9)+ scale_color_manual(values=randpalette14) + 
  theme_classic() + labs(x=paste0("Axis 1 (",Axis1_PctVar, "%)" ), y =paste0("Axis 2 (",Axis2_PctVar, "%)" ) , shape="DFU Healed by 12 Weeks") + coord_equal()

ggsave(PCoa_Plot_By_CC, file="~/Documents/Saureus_Genomics_Paper/PCoa_Accessory_CC.pdf", width=7, height=7)
# Healing group (healed or didn't by 12 weeks) has R^2 of .0348, while CC has R^2 of 0.70269
set.seed(19104)
healinggroup=PCOAresultVectors$HealedBy12
CCgroup=PCOAresultVectors$CCLabel


HealingPermanova = vegan::adonis(formula = JacDist~ healinggroup)
CC_Permanova = vegan::adonis(formula = JacDist~ CCgroup)
both  = vegan::adonis(formula = JacDist~ CCgroup + healinggroup)


#'Misclustered' genomes
#########################
# DORN1352 and DORN1334 are clustered with CC8 instead of the other CC15s
# What are the different genes between them and other CC15s? 
CC15Genes = CCMap %>% filter(CCLabel == "CC15")
CC15Genes$DORN
cc15s = Presence_Absence_Accessory[CC15Genes$DORN, ]
outsiders = cc15s[c("DORN1334", "DORN1352"),]
maincluster = cc15s[setdiff(CC15Genes$DORN,c("DORN1334", "DORN1352")),]

Differential15s = data.frame( OutsiderIncidence= colSums(outsiders), InsiderIncidence = colSums(maincluster))

maincluster = Differential15s %>% filter((OutsiderIncidence==2 & InsiderIncidence==0) | (OutsiderIncidence==0 & InsiderIncidence==10))

OutsiderPresentCC15 = maincluster %>% filter(OutsiderIncidence==2)
OutsiderPresentCC15$Gene = row.names(OutsiderPresentCC15)

OutsiderAbsentcCC15 = maincluster %>% filter(OutsiderIncidence==0)
OutsiderAbsentcCC15$Gene = row.names(OutsiderAbsentcCC15)

outsider_present_phageplasmids = Genes_phage_plasmid %>% filter(X %in% OutsiderPresentCC15$Gene)
outsider_present_phageplasmids$X.1=NULL
outsider_present_phageplasmids$rowsums = rowSums(outsider_present_phageplasmids[2:ncol(outsider_present_phageplasmids)])
outsider_present_phageplasmids %>% filter(rowsums>0) # group_2145, group_600, pepA1, group_1678 seem to be associated with some phages, and group_2219, group_4466, group_505 with others


outsider_absent_phageplasmids = Genes_phage_plasmid %>% filter(X %in% OutsiderAbsentcCC15$Gene)
outsider_absent_phageplasmids$rowsums = rowSums(outsider_absent_phageplasmids[2:ncol(outsider_absent_phageplasmids)])
outsider_absent_phageplasmids %>% filter(rowsums>0) # group_2145, group_600, pepA1, group_1678 seem to be associated with some phages, and group_2219, group_4466, group_505 with others

# DORN1081, DORN1082, DORN1085, DORN1086 are really ST188 which could be considered its own CC 
# What are the different genes between them and other 'CC1s'? 
CC1Genes = CCMap %>% filter(CCLabel == "CC1")
CC1Genes$DORN
cc1s = Presence_Absence_Accessory[CC1Genes$DORN, ]
outsiders = cc1s[c("DORN1081", "DORN1082", "DORN1085", "DORN1086"),]
maincluster = cc1s[setdiff(CC1Genes$DORN,c("DORN1081", "DORN1082", "DORN1085", "DORN1086")),]

Differential1s = data.frame( OutsiderIncidence= colSums(outsiders), InsiderIncidence = colSums(maincluster))

maincluster = Differential1s %>% filter((OutsiderIncidence==4 & InsiderIncidence==0) | (OutsiderIncidence==0 & InsiderIncidence==27))

OutsiderPresent = maincluster %>% filter(OutsiderIncidence==4)
OutsiderAbsent = maincluster %>% filter(OutsiderIncidence==0)

outsider_present_phageplasmids = Genes_phage_plasmid %>% filter(X %in% row.names(OutsiderPresent))
outsider_present_phageplasmids$rowsums = rowSums(outsider_present_phageplasmids[2:ncol(outsider_present_phageplasmids)])
outsider_present_phageplasmids %>% filter(rowsums>0) 

outsider_absent_phageplasmids = Genes_phage_plasmid %>% filter(X %in% row.names(OutsiderAbsent))
outsider_absent_phageplasmids$rowsums = rowSums(outsider_absent_phageplasmids[2:ncol(outsider_absent_phageplasmids)])
outsider_absent_phageplasmids %>% filter(rowsums>0) # group_2145, group_600, pepA1, group_1678 seem to be associated with some phages, and group_2219, group_4466, group_505 with others


# GO annotation of phage and plasmid-associated functions
#########################################################
gene_to_id = RoaryComplete %>% select(Gene, NewpangenomeID)
gene_to_id = gene_to_id %>% left_join(FullAnnotated, by="NewpangenomeID")
Genes_phage_plasmid$Gene = Genes_phage_plasmid$X
AnnotatedPhagePlasmids = Genes_phage_plasmid %>% left_join(gene_to_id, by="Gene")

JustPhages = AnnotatedPhagePlasmids %>% select(c("Gene", "GO.IDs", colnames(AnnotatedPhagePlasmids)[grepl("Phage",colnames(AnnotatedPhagePlasmids))]))

phagenames=colnames(AnnotatedPhagePlasmids)[grepl("Phage",colnames(AnnotatedPhagePlasmids))]
# Mapping genes to GOs
MakeGOMapForCP = AnnotatedPhagePlasmids %>% select("Gene", "GO.IDs")
UncollapsedGOs = tidyr::separate_rows(MakeGOMapForCP,GO.IDs, sep="; " )
UncollapsedGOs = data.frame(UncollapsedGOs)
UncollapsedGOs = UncollapsedGOs %>% filter(GO.IDs!= "" & !is.na(GO.IDs))
termtogene = UncollapsedGOs %>% select(GO.IDs,Gene)
colnames(termtogene) = c( "GOTerms", "Gene")

# Present in phages
present_in_phages = JustPhages %>% select(-GO.IDs)
present_in_phages$total = rowSums(present_in_phages[2:ncol(present_in_phages)])
PhagePresentGenes = (present_in_phages %>% filter(total>0))$Gene
termtogene = termtogene %>% filter(GOTerms!="" & !is.na(GOTerms))
PhageEnriched = enricher(PhagePresentGenes, TERM2GENE = termtogene)

PhageProcessGOsEnriched  = (PhageEnriched@result %>% filter(p.adjust<.05) %>% select(ID, p.adjust,geneID) %>% arrange(ID,p.adjust ))

PhageProcessGOsEnriched$ID

PhageProcessGOsEnriched$IDNoOntology = sapply(PhageProcessGOsEnriched$ID, function(x) paste(str_split(string=x, pattern=":")[[1]][2:3], collapse=":"))

PhageProcessGOsEnriched$Description = sapply(PhageProcessGOsEnriched$IDNoOntology, function(x) toString(go2term(x)["Term"]))

PhageProcessGOsEnriched$geneID = sapply(PhageProcessGOsEnriched$geneID, function(x) str_replace_all(x,"/",";"))

PhageProcessGOsEnriched = PhageProcessGOsEnriched %>% select(ID, IDNoOntology, Description, p.adjust,geneID)
write.csv(PhageProcessGOsEnriched, file="~/Documents/DataInputGithub/data/IntraPatient/Phages/EnrichedGOsPhageGenes.csv")

# Present in plasmids
JustPlasmids = AnnotatedPhagePlasmids %>% select(-phagenames)
Present_In_Plasmids = JustPlasmids %>% select(-c(X, NewpangenomeID, GO.IDs))
Present_In_Plasmids$rowsums = rowSums(Present_In_Plasmids[1:(ncol(Present_In_Plasmids)-2)])
PlasmidPresentGenes = (Present_In_Plasmids %>% filter(rowsums > 0))$Gene

PlasmidEnriched = enricher(PlasmidPresentGenes, TERM2GENE = termtogene)

PlasmidEnrichedGOs  = (PlasmidEnriched@result %>% filter(p.adjust<.05) %>% select(ID, p.adjust,geneID) %>% arrange(ID,p.adjust ))

PlasmidEnrichedGOs$IDNoOntology = sapply(PlasmidEnrichedGOs$ID, function(x) paste(str_split(string=x, pattern=":")[[1]][2:3], collapse=":"))

PlasmidEnrichedGOs$Description = sapply(PlasmidEnrichedGOs$IDNoOntology, function(x) toString(go2term(x)["Term"]))

PlasmidEnrichedGOs$geneID = sapply(PlasmidEnrichedGOs$geneID, function(x) str_replace_all(x,"/",";"))

PlasmidEnrichedGOsWrite = PlasmidEnrichedGOs %>% select(ID, IDNoOntology, Description, p.adjust,geneID)


PlasmidEnrichedGOs %>% filter(IDNoOntology %in% intersect(PlasmidEnrichedGOs$IDNoOntology, PhageProcessGOsEnriched$IDNoOntology))


PlasmidEnrichedGOs %>% filter(!(IDNoOntology %in% intersect(PlasmidEnrichedGOs$IDNoOntology, PhageProcessGOsEnriched$IDNoOntology)))
(PhageProcessGOsEnriched %>% filter(!(IDNoOntology %in% intersect(PlasmidEnrichedGOs$IDNoOntology, PhageProcessGOsEnriched$IDNoOntology))))$Description

dim(PlasmidEnrichedGOsWrite)

write.csv(PlasmidEnrichedGOsWrite, "~/Documents/DataInputGithub/data/IntraPatient/Plasmids/EnrichedGOsPlasmidGenes.csv")

PhagesNotPlasmids = setdiff(PhageProcessGOsEnriched$ID, PlasmidEnrichedGOs$ID)

PlasmidsNotPhages = setdiff(PlasmidEnrichedGOs$ID,PhageProcessGOsEnriched$ID)


RoaryOutputPresenceAbsence = RoaryOutputPresenceAbsence %>% left_join(FullAnnotated,by="NewpangenomeID")
RoaryOutputPresenceAbsence = RoaryOutputPresenceAbsence %>% select(-OldPangenomeID)

DORNList = colnames(RoaryOutputPresenceAbsence)[which(grepl("DORN", colnames(RoaryOutputPresenceAbsence)))]
patient_genome = patient_genome %>% filter(DORN %in% DORNList)

# Presence/absence of each gene by patient 
GenePresenceDB = RoaryOutputPresenceAbsence %>% select(Gene,NewpangenomeID, GO.IDs)

pat="141"


for(pat in unique(patient_genome$patient_id)){
  GenomesPatient= (patient_genome %>% filter(patient_id==pat))$DORN
  ColsPatient = RoaryOutputPresenceAbsence %>% dplyr::select(NewpangenomeID,GO.IDs, GenomesPatient)
  if(length(GenomesPatient) == 1){
    ColsPatient$Sums = ColsPatient[,3]
    
  }else{
    TestFrame=ColsPatient[,3:ncol(ColsPatient)]#apply(ColsPatient[,3:ncol(ColsPatient)],2, function(x) as.numeric(x))
    ColsPatient$Sums = rowSums(TestFrame)
    
  }
  
  ColsPatient$InPatient = if_else(ColsPatient$Sums>0, 1,0)
  GenePresenceDB[paste0("patient_", pat)] = ColsPatient$InPatient
}
write.csv(GenePresenceDB, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient.csv")

# Gene by patient (CC5 genomes only)
#####################################
CC5IDs = (CCMap %>% filter(CCLabel=="CC5"))$DORN

patient_genome_CC5 = patient_genome %>% filter(DORN %in% CC5IDs)
GenePresenceDB_CC5 = RoaryOutputPresenceAbsence %>% select(Gene,NewpangenomeID, GO.IDs)

for(pat in unique(patient_genome_CC5$patient_id)){
  GenomesPatient= (patient_genome_CC5 %>% filter(patient_id==pat))$DORN
  ColsPatient = RoaryOutputPresenceAbsence %>% dplyr::select(NewpangenomeID,GO.IDs, GenomesPatient)
  if(length(GenomesPatient) == 1){
    ColsPatient$Sums = ColsPatient[,3]
    
  }else{
    TestFrame=ColsPatient[,3:ncol(ColsPatient)]#apply(ColsPatient[,3:ncol(ColsPatient)],2, function(x) as.numeric(x))
    ColsPatient$Sums = rowSums(TestFrame)
    
  }
  
  ColsPatient$InPatient = if_else(ColsPatient$Sums>0, 1,0)
  GenePresenceDB_CC5[paste0("patient_", pat)] = ColsPatient$InPatient
}
write.csv(GenePresenceDB_CC5, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient_CC5.csv")




# Gene by patient (CC1 genomes only)
#####################################
CC1IDs = (CCMap %>% filter(CCLabel=="CC1"))$DORN

patient_genome_CC1 = patient_genome %>% filter(DORN %in% CC1IDs)
GenePresenceDB_CC1 = RoaryOutputPresenceAbsence %>% select(Gene,NewpangenomeID, GO.IDs)

for(pat in unique(patient_genome_CC1$patient_id)){
  GenomesPatient= (patient_genome_CC1 %>% filter(patient_id==pat))$DORN
  ColsPatient = RoaryOutputPresenceAbsence %>% dplyr::select(NewpangenomeID,GO.IDs, GenomesPatient)
  if(length(GenomesPatient) == 1){
    ColsPatient$Sums = ColsPatient[,3]
    
  }else{
    TestFrame=ColsPatient[,3:ncol(ColsPatient)]#apply(ColsPatient[,3:ncol(ColsPatient)],2, function(x) as.numeric(x))
    ColsPatient$Sums = rowSums(TestFrame)
    
  }
  
  ColsPatient$InPatient = if_else(ColsPatient$Sums>0, 1,0)
  GenePresenceDB_CC1[paste0("patient_", pat)] = ColsPatient$InPatient
}
write.csv(GenePresenceDB_CC1, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient_CC1.csv")



# Gene by patient (CC8 genomes only)
#####################################
CC8IDs = (CCMap %>% filter(CCLabel=="CC8"))$DORN

patient_genome_CC8 = patient_genome %>% filter(DORN %in% CC8IDs)
GenePresenceDB_CC8 = RoaryOutputPresenceAbsence %>% select(Gene,NewpangenomeID, GO.IDs)

for(pat in unique(patient_genome_CC8$patient_id)){
  GenomesPatient= (patient_genome_CC8 %>% filter(patient_id==pat))$DORN
  ColsPatient = RoaryOutputPresenceAbsence %>% dplyr::select(NewpangenomeID,GO.IDs, GenomesPatient)
  if(length(GenomesPatient) == 1){
    ColsPatient$Sums = ColsPatient[,3]
    
  }else{
    TestFrame=ColsPatient[,3:ncol(ColsPatient)]#apply(ColsPatient[,3:ncol(ColsPatient)],2, function(x) as.numeric(x))
    ColsPatient$Sums = rowSums(TestFrame)
    
  }
  
  ColsPatient$InPatient = if_else(ColsPatient$Sums>0, 1,0)
  GenePresenceDB_CC8[paste0("patient_", pat)] = ColsPatient$InPatient
}
write.csv(GenePresenceDB_CC8, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GeneByPatient_CC8.csv")



Accessory_Genes_by_Patient = GenePresenceDB 
Accessory_Genes_by_Patient$Sums = rowSums(Accessory_Genes_by_Patient[, 4:ncol(Accessory_Genes_by_Patient)])
Accessory_Genes_by_Patient = Accessory_Genes_by_Patient %>% filter(Sums<60)

# Core genes by patient (genes all patients had at least one S aureus genome with)
CoreGenes = RoaryOutputPresenceAbsence %>%  filter(NewpangenomeID %in% setdiff(GenePresenceDB$NewpangenomeID, Accessory_Genes_by_Patient$NewpangenomeID))


FullListGOs = RoaryOutputPresenceAbsence$GO.IDs
FullListGOs = (sapply(FullListGOs, function(x) str_split(x, pattern="; ")  ))
allGOs = unique(Reduce(c, FullListGOs))
Functions_Processes_Only = allGOs[!grepl("C:",allGOs)]
Functions_Processes_Only=Functions_Processes_Only[Functions_Processes_Only!=""]
Functions_Processes_Only = Functions_Processes_Only[!is.na(Functions_Processes_Only)]

Accessory_Genes_by_Patient$Sums=NULL
PlotDF = data.frame()

Functions_Processes_Only = Functions_Processes_Only[!is.na(Functions_Processes_Only)]

for(p in unique(patient_genome$patient_id)){
  patientrow=c(p)
  print(p)
  for(g in Functions_Processes_Only){
   
    PatientDB= Accessory_Genes_by_Patient[Accessory_Genes_by_Patient[,paste0("patient_", p)]==1 , ]
    Occurences = length(grep(g, PatientDB$GO.IDs))
    patientrow = append(patientrow, Occurences)
  }
  PlotDF = rbind(PlotDF, patientrow)
}


colnames(PlotDF) = c("PatientID", Functions_Processes_Only)

write.csv(PlotDF, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GOTermOccurences_By_Patient.csv")


PlotDF = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/GOTermOccurences_By_Patient.csv")
AllTerms = colnames(PlotDF)[!(colnames(PlotDF) %in% c("X", "PatientID"))]

write.table(AllTerms, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/FullGOaccessorylist.txt", quote=F, col.names=F, row.names=F)

# Add to the big plot 
######################
PatientInfo = read.csv("~/Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")

PatientInfo$DORN=paste0("DORN", PatientInfo$Doern.lab.bank.)
PatientInfo = PatientInfo %>% select(DORN, patient_id, visit,)
IncludedIsolates = CCMap %>% filter(!(DORN %in% c("DORN429","DORN1176", "DORN685", "DORN946")))
PatientsCC = IncludedIsolates %>% left_join(PatientInfo, by="DORN")
PatientsCC$Week = PatientsCC$visit*2 
PatientsCC$amount=1
ByCC = PatientsCC %>% group_by(patient_id, CCLabel) %>% summarize(Total=sum(amount)) %>% unique()

listBigTree= c("CC1", "CC12","CC133", "CC15", "CC20", "CC22", "CC30", "CC398", "CC45", "CC5", "CC59", "CC7", "CC72", "CC8", "CC80", "CC9", "CC97" )

Indices_Remove = which(listBigTree %in% setdiff(listBigTree, ByCC$CCLabel))


randpalette14

WeekHealed$patient_id = WeekHealed$patient
ByCCweek = ByCC %>% left_join(WeekHealed %>% select(patient_id, week_healed)) %>% unique()
ByCCweek = ByCCweek %>% mutate(WeekNumeric = if_else(is.na(week_healed), 100, as.integer(week_healed)))

CCplot = ggplot(ByCCweek, aes(x=factor(patient_id),y=Total, fill=CCLabel)) + geom_bar(stat="identity") + scale_fill_manual(values=randpalette14) + theme_classic() + theme(axis.text.x=element_text(angle=90))

CCplot$data$patient_id = factor(CCplot$data$patient_id , levels= unique((ByCCweek %>% arrange(week_healed))$patient_id))

OrderPlots = unique((ByCCweek %>% arrange(week_healed))$patient_id)

RoaryComplete




# Make accessory gene funciton plot 
# Only 491 GO terms are actually found in 
Sums = colSums(PlotDF[2:ncol(PlotDF)])
Sums = Sums[Sums>0]

Included = PlotDF %>% select(PatientID, names(Sums))

IncludedTerms=colnames(Included)[colnames(Included)!="PatientID"]
write.table(IncludedTerms, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/IncludedGOaccessorylist.txt", quote=F, col.names=F, row.names=F)
GOterms = read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/GOslimTermsAccessory.txt", sep="\t", header=F)
GOterms = GOterms %>% select(V1,V2,V3)

GOtermsList = (sapply(GOterms$V3, function(x) str_split(x, pattern=";")  ))

Represented = unique(Reduce(c, GOtermsList))

SlimList=c()
goList = c()
for(GO in Represented){
  print(GO)
  Selected=FALSE
  
  # Start at most specific/rarest GOslim term
  TestRow = nrow(GOterms)
  while((Selected==FALSE) & (TestRow>0)){
    if(grepl(pattern=GO,x=GOterms[TestRow, "V3"])){
      #print(paste0("Assign ",  GO, " to ",GOterms[TestRow, "V1"]))
      goList = append(goList, GO)
      SlimList = append(SlimList,GOterms[TestRow, "V1"] )
      Selected=TRUE
    }
    # go a little less specific
    TestRow=TestRow-1
  }
  
}





#"GO:0008199"


MapGOstoSlim = data.frame(SlimTerm = SlimList, GoTerm=goList)
MyGOlist = data.frame(MyGO = IncludedTerms)
MyGOlist$GoTerm = sapply(MyGOlist$MyGO, function(x) str_replace(x,"P.", ""))
MyGOlist$GoTerm = sapply(MyGOlist$GoTerm, function(x) str_replace(x,"F.", ""))
MyGOlist$GoTerm = sapply(MyGOlist$GoTerm, function(x) str_replace(x,"\\.", ":"))

MyGOlist = MyGOlist %>% left_join(MapGOstoSlim, by="GoTerm")


OrderGOs = (MyGOlist %>% arrange(SlimTerm))$MyGO

IncludedPresAbsence = Included[,2:ncol(Included)]
IncludedPresAbsence[IncludedPresAbsence>0] <- 1
IncludedPresAbsence$PatientID = Included$PatientID

write.table(Included, )
MeltedGO = IncludedPresAbsence %>% reshape2::melt(id.vars="PatientID")
MeltedGO$Type = if_else(grepl(pattern="P:",MeltedGO$variable), "Process", "Function")

GPlot = ggplot(MeltedGO, aes(x=factor(PatientID), y=variable, fill=factor(value))) + geom_tile() +scale_fill_manual(values=c("white","aquamarine3")) 

GPlot$data$PatientID = factor(GPlot$data$PatientID, levels=OrderPlots)
GPlot$data$variable = factor(GPlot$data$variable, levels=sort(unique(GPlot$data$variable)))

GPlotG = GPlot + theme(axis.text.x=element_text(angle=90))

GPlotG$data$variable = factor(GPlotG$data$variable, levels=OrderGOs)
ggsave(GPlotG, height=30, width=15, file="Documents/Saureus_Genomics_Paper/AccessoryGOs.pdf")


GOterms$SlimTerm = GOterms$V1
MyGOlist = MyGOlist %>% left_join(GOterms, by="SlimTerm") 
MyGOlist = MyGOlist %>% select(V2, GoTerm, SlimTerm) %>% arrange(SlimTerm)
write.csv(MyGOlist, file="/Users/amycampbell/Documents/DataInputGithub/data/SlimTermFunctions.csv")



HealedPatientsGO = IncludedPresAbsence %>% filter(PatientID %in% HealedPatients$patient) 
UnhealedPatientsGO = IncludedPresAbsence %>% filter(PatientID %in% UnhealedPatients$patient) 
HealedPatientsGO$PatientID = NULL
UnhealedPatientsGO$PatientID = NULL
PresentHealed = names(colSums(HealedPatientsGO)[colSums(HealedPatientsGO)>0])
PresentUnHealed = names(colSums(UnhealedPatientsGO)[colSums(UnhealedPatientsGO)>0])

UnhealedExclusive= setdiff(PresentUnHealed, PresentHealed)

sapply(UnhealedExclusive, function(x) str_replace_all(x, "", ":"))

Accessory_Genes_by_Patient
UnhealedExcl = sapply(UnhealedExclusive, function(x) str_replace_all(x, "\\.", ":"))
UnhealedExcl = sapply(UnhealedExcl, function(x) str_replace_all(x, "F:", ""))
UnhealedExcl = sapply(UnhealedExcl, function(x) str_replace_all(x, "P:", ""))
for (u in UnhealedExcl){
  print(u)
}
  print(UnhealedExcl)
# 
# Blast2GoAnnotationFile$GeneName = GeneNames
# 
# 
# Blast2GoAnnotationFileNoComponents = Blast2GoAnnotationFile %>% filter( grepl(GO.IDs, pattern="F:") | grepl(GO.IDs, pattern="P:")) %>%filter( !is.na(X.GO))
# Blast2GoAnnotationFileNoComponentsAccessory = Blast2GoAnnotationFileNoComponents %>% filter(GeneName  %in% AccessoryOutput$Gene)
# 
# DORN925Genes = (RoaryOutputPresenceAbsence %>% filter(DORN925==1))$Gene
# DORN925genes_Annotations = Blast2GoAnnotationFile %>% filter (GeneName %in% DORN925Genes)
# 
# dim(DORN925genes_Annotations %>% filter( grepl(GO.IDs, pattern="F:") | grepl(GO.IDs, pattern="P:")) %>%filter( !is.na(X.GO)))
# dim(DORN925genes_Annotations  %>%filter( !is.na(X.GO)))
# 
# # equivalent to not using unique() since theres one line per sequence name 
# length(unique((Blast2GoAnnotationFile %>% filter(!is.na(X.GO)))$SeqName)) # 3768
# length(unique((Blast2GoAnnotationFile %>% filter((Enzyme.Codes !="")))$SeqName)) # 208 
# length(unique((Blast2GoAnnotationFile %>% filter(!(InterPro.GO.IDs %in% c("no IPS match", "no GO terms"))))$SeqName)) # 3003 
# # GO is definitely the move, therefore. 
# 
# Blast2GoAnnotationFileNoComponents = Blast2GoAnnotationFile %>% filter( grepl(GO.IDs, pattern="F:") | grepl(GO.IDs, pattern="P:"))
# length(unique((Blast2GoAnnotationFileNoComponents %>% filter(!is.na(X.GO)))$SeqName))







# 
# GOIDs_Expanded = Blast2GoAnnotationFile$GO.IDs
# length(GOIDs_Expanded)
# allGOinstances=c()
# for(go in 1:length(GOIDs_Expanded)){
#   golist = GOIDs_Expanded[[go]]
#   print(length((stringr::str_split(golist, pattern="; "))[[1]]))
#   allGOinstances = append(allGOinstances,(stringr::str_split(golist, pattern="; "))[[1]])
# }
# 
# allGOinstances = allGOinstances[allGOinstances!=""]
# GOtable = table(allGOinstances)
# GOterm = names(GOtable)
# GOdf = data.frame(GOtable)
# colnames(GOdf) = c("Term", "Count")
# GOdf = GOdf %>% arrange(-Count)d
# GOdf = GOdf %>% mutate(GO_Type = case_when(grepl(Term, pattern="C:") ~ "Cellular Component", 
#                                             grepl(Term, pattern="P:") ~ "Biological Process",
#                                             grepl(Term, pattern="F:") ~ "Molecular Function"))
# #GOdfNoComponent = GOdf %>% filter(!grepl(Term, pattern="C:"))
# GOdfAtLeast5 = GOdf %>% filter(Count >10)
# BarPlotGOs = ggplot(GOdfAtLeast5, aes(x=Term, y=Count, fill=GO_Type)) + geom_bar( stat="identity")
# BarPlotGOs$data$Term = factor(BarPlotGOs$data$Term, levels=GOdfAtLeast5$Term) 
# BarPlotGOs = BarPlotGOs + scale_fill_manual(values=c("#B53737","#190BDA","#E2D629" )) + theme_classic() + theme(axis.text.x=element_text(size=6.5, angle=70, hjust=.4, vjust=.4), plot.title=element_text(size=15, hjust=.5, face="bold")) + ggtitle("GO Terms with >10 Frequency in the Pan-Genome")
# #ggsave(BarPlotGOs, width=12,height=6, file="data/Phylogeny2022Data/GOTermsAtLeast10.pdf")
# 
# 
# GOdfAtLeast10_Fxns  = GOdfAtLeast5 %>% filter(GO_Type != "Cellular Component")
# BarPlotGOsFxns = ggplot(GOdfAtLeast10_Fxns, aes(x=Term, y=Count, fill=GO_Type)) + geom_bar( stat="identity",color="black") + theme_classic() + theme(axis.text.x=element_text(size=8, angle=70, hjust=.4, vjust=.4), plot.title=element_text(size=15, hjust=.5, face="bold")) + scale_fill_manual(values=c("#190BDA","#E2D629" )) + ggtitle("Functional GO Terms with >10 Frequency in the Pan-Genome")
# BarPlotGOsFxns$data$Term = factor(BarPlotGOsFxns$data$Term, levels = (GOdfAtLeast10_Fxns %>% arrange(-Count))$Term)
# 
# ggsave(BarPlotGOsFxns, width=13,height=7, file="data/Phylogeny2022Data/GOTermsAtLeast10FunctionsOnly.pdf")

write.csv(names(Included)[names(Included)!="PatientID"], )


IncludedPresAbsence118 = IncludedPresAbsence %>% filter(PatientID == "118")
IncludedPresAbsence118$PatientID = NULL
IncludedPresAbsenceOthers = IncludedPresAbsence %>% filter(PatientID !="118")
IncludedPresAbsenceOthers$PatientID = NULL



GenePatient = Accessory_Genes_by_Patient$Gene
PatientBYGene = t(Accessory_Genes_by_Patient[, 4:ncol(Accessory_Genes_by_Patient)])
patientNames = rownames(PatientBYGene)

PatientBYGene = data.frame(PatientBYGene)

colnames(PatientBYGene) = GenePatient

PatientBYGene$patient = sapply(patientNames, function(x) as.numeric(as.character(str_split(x, "patient_")[[1]][2])))

PatientBYGene = PatientBYGene %>% left_join(WeekHealed %>% select(HealedBy12,patient ) , by='patient')
PatientBYGene$patient=NULL

genes = colnames(PatientBYGene)[colnames(PatientBYGene)!="HealedBy12"]

GeneSums = colSums(PatientBYGene[,genes])
GeneSums = GeneSums[GeneSums>0]
genes = names(GeneSums)
pvals = c()
OR=c()
for(g in genes){
  print(g)
  FishTest=fisher.test(table(PatientBYGene[,c(g, "HealedBy12")]))
  
  pvals=append(pvals, FishTest$p.value)
  OR=append(OR, FishTest$estimate)
  
}

FullDF = data.frame(Pval=pvals, OddsRatio=OR, gene=genes)
FullDF$Padj = p.adjust(FullDF$Pval, method="BH")

HigherInHealed = FullDF %>% filter(Padj <.05 & OddsRatio>1)
HigherInUnHealed = FullDF %>% filter(Padj <.05 & OddsRatio<1)

# when higher in 'yes healed' than 'no healed', OR>1 
# when higher in 'no healed' than 'yes healed', OR<1 
MoreCommon_Unhealed_Genes = HigherInUnHealed$gene
MoreCommon_Healed_Genes = HigherInHealed$gene

UnhealedEnriched = Accessory_Genes_by_Patient %>% filter(Gene %in% MoreCommon_Unhealed_Genes)
HealedEnriched = Accessory_Genes_by_Patient %>% filter(Gene %in% MoreCommon_Healed_Genes)


GosUnhealed = UnhealedEnriched %>% select(GO.IDs,Gene) %>% filter(!is.na(GO.IDs) & !(GO.IDs==""))
GosUnhealed$GO.IDs = sapply(GosUnhealed$GO.IDs, function(x) str_remove_all(x," "))
GosUnhealed = data.frame(tidyr::separate_rows(GosUnhealed, GO.IDs, sep=";"))
colnames(GosUnhealed) = c("GO_terms","Gene")


GosHealed = HealedEnriched %>% select(GO.IDs,Gene) %>% filter(!is.na(GO.IDs) & !(GO.IDs==""))
GosHealed$GO.IDs = sapply(GosHealed$GO.IDs, function(x) str_remove_all(x," "))
GosHealed = data.frame(tidyr::separate_rows(GosHealed, GO.IDs, sep=";"))
colnames(GosHealed) = c("GO_terms","Gene")

resultsHealed= clusterProfiler::enricher(GosHealed$Gene, TERM2GENE = GosHealed)
resultsUnhealed = clusterProfiler::enricher(UnhealedEnriched$Gene, TERM2GENE = GosUnhealed)


GoMap = Accessory_Genes_by_Patient %>% select(GO.IDs,Gene) %>% filter(!is.na(GO.IDs) & !(GO.IDs==""))
GoMap$GO.IDs = sapply(GoMap$GO.IDs, function(x) str_remove_all(x," "))
GoMap = data.frame(tidyr::separate_rows(GoMap, GO.IDs, sep=";"))
colnames(GoMap) = c("GO_terms","Gene")

GSEAdf = FullDF$OddsRatio
GSEAdf = log2(GSEAdf)
names(GSEAdf)= FullDF$gene
GSEAdf = GSEAdf[!is.infinite(GSEAdf)]
ResultsGSEA= clusterProfiler::GSEA(rev(sort(GSEAdf)), TERM2GENE = GoMap)
resultsUnhealed = clusterProfiler::GSEA(rev(sort(UnhealedEnriched$Gene)), TERM2GENE = GosUnhealed)



