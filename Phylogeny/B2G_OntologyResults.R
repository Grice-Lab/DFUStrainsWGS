library(dplyr)
library(stringr)
library(ggplot2)
#Blast2GoMappingFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoMappingResults_Original.txt", sep="\t",header=T,row.names=NULL)
#Blast2GoMappingFileNew = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoMappingResults_2149.txt", sep="\t",header=T,row.names=NULL)
randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#F6BE00",
                "#33A02C","#E6F5C9")

Patient_Genome_Info = read.csv("~/Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")
Patient_Genome_Info
Patient_Genome_Info$DORN = paste0("DORN", Patient_Genome_Info$Doern.lab.bank.)
patient_genome = Patient_Genome_Info %>% select(patient_id, DORN)
WeekHealed = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/staphyloxanthin_paper_data.csv") %>% select(patient, week_healed)

CCMap = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Phylogeny2022Data/CCMapPlotting.csv")

Blast2GoAnnotationFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoAnnotations_Original.txt", sep="\t",header=T,row.names=NULL) %>% select(GO.IDs, SeqName)
Blast2GoMapFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoMappingResults_Original.txt", sep="\t",header=T,row.names=NULL)# %>% select(GO.IDs, SeqName)

Blast2GoAnnotationFileNew = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoAnnotations_2149.txt", sep="\t",header=T,row.names=NULL) %>% select(GO.IDs, SeqName)

Gos =c("GO:0006457",
"GO:0140662",
"GO:0046677",
"GO:0019835",
"GO:0140664",
"GO:0140359",
"GO:0016987",
"GO:0051301",
"GO:0051715",
"GO:0007155",
"GO:0009372",
"GO:0140663",
"GO:0031647",
"GO:0043462",
"GO:0043937",
"GO:0042026",
"GO:0019836",
"GO:0043708",
"GO:0009294",
"GO:0090729",
"GO:0032196",
"GO:0046685",
"GO:0140658",
"GO:0106256",
"GO:0006276")
Blast2GoAnnotationFile %>% filter()
# SeqNames for the 'new' run of blast2go are the same as seqnames in most recent roary run
Blast2GoAnnotationFileNew$NewpangenomeID = Blast2GoAnnotationFileNew$SeqName
# But seqnames for the 'old' run of blast2go are the old roary run
Blast2GoAnnotationFile$OldPangenomeID = Blast2GoAnnotationFile$SeqName

OldToNewMap=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Map_New_To_Old_Pangenome.csv")

Blast2GoAnnotationFileNew = Blast2GoAnnotationFileNew %>% left_join(OldToNewMap, by="NewpangenomeID")
Blast2GoAnnotationFile = Blast2GoAnnotationFile %>% left_join(OldToNewMap, by="OldPangenomeID") %>% filter(!is.na(NewpangenomeID))

Blast2GoAnnotationFileNew = Blast2GoAnnotationFileNew %>% select(GO.IDs,NewpangenomeID)
Blast2GoAnnotationFile = Blast2GoAnnotationFile %>% select(GO.IDs, NewpangenomeID)
FullAnnotated  = rbind(Blast2GoAnnotationFileNew, Blast2GoAnnotationFile)

#Blast2GoMappingFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoAnnotations_Original.txt", sep="\t",header=T,row.names=NULL)
RoaryOutput = read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new.csv", sep=',')
RoaryOutput$NewpangenomeID = 1:nrow(RoaryOutput)

# For each item in OldToNewMap$NewpangenomeID, I want to figure out which row has an entry that == it and then stick it into a variable in that row
#for(id in OldToNewMap$NewpangenomeID){
  #listobj=apply(RoaryOutput, 2, function(x) which(x == id))
 # Index = (which(apply(RoaryOutput, 1, function(x) any(x==id))))
#  RoaryOutput[Index,"NewpangenomeID"] = id
#}
#RoaryOutput = RoaryOutput %>% left_join(OldToNewMap, by="NewpangenomeID")
#write.csv(RoaryOutput, file="/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv")


RoaryComplete = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv")
RoaryComplete
RoaryOutputPresenceAbsence = RoaryComplete
RoaryOutputPresenceAbsence = RoaryOutputPresenceAbsence %>% select(Gene, OldPangenomeID, NewpangenomeID, colnames(RoaryOutputPresenceAbsence)[which(grepl("DORN", colnames(RoaryOutputPresenceAbsence)))])


RoaryOutputPresenceAbsence = RoaryOutputPresenceAbsence %>% mutate_at(vars(contains("DORN")), function(x) if_else(x=="", 0,1))

#RoaryOutputPresenceAbsence[,4:ncol(RoaryOutputPresenceAbsence)] = apply(RoaryOutputPresenceAbsence[,4:ncol(RoaryOutputPresenceAbsence)], 1, function(x) as.numeric(as.character(x)))


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

randpalette14 = randpalette18[-c(3,6,8,15)]

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

WeekHealed$week_healed
WeekHealed$HealedBy12 = if_else(WeekHealed$week_healed>12 | is.na(WeekHealed$week_healed), "No", "Yes")
HealedPatients = WeekHealed %>% filter(HealedBy12=="Yes")
UnhealedPatients = WeekHealed %>% filter(HealedBy12=="No")


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



