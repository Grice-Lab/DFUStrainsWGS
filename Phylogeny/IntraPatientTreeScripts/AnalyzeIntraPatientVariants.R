library(dplyr)
library(stringr)
library(ggplot2)
# Amy campbell
# Analyzing the intra-patient gene variant contents 

# CC color scheme
#################
cc_colors = read.csv("~/Documents/DataInputGithub/data/Phylogeny2022Data/CC_Colorscheme.csv")

# Gene ortholog presence/absence from Roary
##############################################
GeneAnnotations = read.csv("~/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv") %>% select(Gene,Annotation)

# Gene ortholog presence/absence by phage and plasmid
###########################################################
genepresence_phage_plasmid = read.csv("~/Documents/DataInputGithub/data/IntraPatient/GenePresence_ByPhage_Plasmid_Updated_Hclust.csv")
genepresence_phage_plasmid$X.1 =NULL

# Presence of the 13 different phage clusters by genome
###########################################################
phagepresence_isolates = read.csv("~/Documents/DataInputGithub/data/IntraPatient/Phages/PresenceAbsenceHClustPhages.csv")

# SNPs / indels that were found to be consistent in any comparisons
##################################################################
SNPs = read.csv("~/Documents/DataInputGithub/data/IntraPatient/SNP_Variants_Consistent.csv")

SNPsBackup = SNPs
CNVs = read.csv("~/Documents/DataInputGithub/data/IntraPatient/CNV_Variants_Consistent.csv")
GeneAndHGT = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/IntraPatientGeneticComparisons/"
colorScheme = read.csv("~/Documents/DataInputGithub/data/IntraPatient/SpectralColorsVariants.csv")# %>% select(TypePresenceCombo,HexVal )

SNPs = SNPs %>% select(Gene, Phenotype, CC, Patient, Cluster1, Cluster2, Type, PresentOrAbsentHigh)
CNVs = CNVs %>% select(Gene, Phenotype, CC, Patient, Cluster1, Cluster2, Type, PresentOrAbsentHigh)
colnames(CNVs) = c("VariantLocus", "Phenotype", "CC", "Patient","Cluster1", "Cluster2", "VariantType", "PresentOrAbsentHigh")
colnames(SNPs) = c("VariantLocus", "Phenotype", "CC", "Patient","Cluster1", "Cluster2", "VariantType", "PresentOrAbsentHigh")

AllVariantsDF = rbind(CNVs, SNPs)


for(filename in list.files(GeneAndHGT)){
  if(file.size(paste0(GeneAndHGT, filename)) > 3){
    comparison = read.csv(paste0(GeneAndHGT, filename))
    comparison = comparison %>% mutate(Cluster1 = min(HighCluster, LowCluster))
    comparison = comparison %>% mutate(Cluster2 = max(HighCluster, LowCluster))
    
    comparison = comparison %>% select(VariantID, phenotype, CCgroup, Patient, Cluster1, Cluster2, VariantType,PresentOrAbsent_in_High)
    comparison$PresentOrAbsent_in_High = sapply(comparison$PresentOrAbsent_in_High, function(x) tolower(x))
    colnames(comparison) = c("VariantLocus", "Phenotype", "CC", "Patient","Cluster1", "Cluster2", "VariantType", "PresentOrAbsentHigh")
    
    AllVariantsDF = rbind(AllVariantsDF, comparison)
    } 
  
}
AllVariantsDF = AllVariantsDF %>% mutate(Comparison = paste(Phenotype, CC, Patient, Cluster1, Cluster2, sep="_"))

AllVariantsDF = AllVariantsDF %>% mutate(Comparison_Gene = paste(Comparison, VariantLocus, sep="_"))
NumberVars = AllVariantsDF %>% filter(!(VariantType %in% c("Synonymous","Phage", "Plasmid"))) %>% group_by(Comparison) %>% summarize(NumVariants = n())
NumberSequenceVars = AllVariantsDF %>% filter(!(VariantType %in% c("Phage", "Plasmid","Gene", "Synonymous"))) %>% group_by(Comparison) %>% summarize(NumSequenceVariants = n())

AllVariantsDF = AllVariantsDF %>% left_join(NumberVars, by="Comparison")
AllVariantsDF = AllVariantsDF %>% left_join(NumberSequenceVars, by="Comparison")
AllVariantsDF$NumSequenceVariants[is.na(AllVariantsDF$NumSequenceVariants)] <- 0
AllVariants_HGTs = AllVariantsDF %>% filter(VariantType %in% c("Phage", "Plasmid"))
AllVariants_SingleGene = AllVariantsDF %>% filter(!(VariantType %in% c("Phage", "Plasmid")))
AllVariants_NonGene = AllVariantsDF %>% filter(VariantType!="Gene")


# 
# MostCommonByGene =
# MostCommonByGene$FixedGeneNames = sapply(MostCommonByGene$VariantLocus, function(x) str_replace(x,"\\'", "_"))
# MostCommonByGene$FixedGeneNames = sapply(MostCommonByGene$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
# MostCommonByGene$FixedGeneNames = sapply(MostCommonByGene$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
# MostCommonByGene$FixedGeneNames = sapply(MostCommonByGene$FixedGeneNames, function(x) str_replace(x,":", "_"))
# MostCommonByGene$FixedGeneNames = sapply(MostCommonByGene$FixedGeneNames, function(x) str_replace(x,"-", "_"))
# MostCommonByGene$Gene = MostCommonByGene$FixedGeneNames 



genepresence_phage_plasmid$FixedGeneNames = sapply(genepresence_phage_plasmid$X, function(x) str_replace(x,"\\'", "_"))
genepresence_phage_plasmid$FixedGeneNames = sapply(genepresence_phage_plasmid$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
genepresence_phage_plasmid$FixedGeneNames = sapply(genepresence_phage_plasmid$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
genepresence_phage_plasmid$FixedGeneNames = sapply(genepresence_phage_plasmid$FixedGeneNames, function(x) str_replace(x,":", "_"))
genepresence_phage_plasmid$FixedGeneNames = sapply(genepresence_phage_plasmid$FixedGeneNames, function(x) str_replace(x,"-", "_"))

# Gene presence by phage
genepresence_phage = genepresence_phage_plasmid %>% select(FixedGeneNames, colnames(genepresence_phage_plasmid)[grepl(colnames(genepresence_phage_plasmid), pattern="Phage")])
genepresence_phage$NumPhages=rowSums(genepresence_phage[2:ncol(genepresence_phage)])
Present_In_Phage = genepresence_phage %>% filter(NumPhages>0) %>% select(FixedGeneNames, NumPhages)

# Gene presence by plasmid
genepresence_plasmids= genepresence_phage_plasmid %>% select(FixedGeneNames, setdiff(colnames(genepresence_phage_plasmid), colnames(genepresence_phage)))
genepresence_plasmids$X=NULL
genepresence_plasmids$NumPlasmids = rowSums(genepresence_plasmids[2:(ncol(genepresence_plasmids)-2)])
Present_In_Plasmids = genepresence_plasmids %>% filter(NumPlasmids>0) %>% select(FixedGeneNames, NumPlasmids)


GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$Gene, function(x) str_replace(x,"\\'", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,":", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,"-", "_"))


test = AllVariants_SingleGene %>% mutate(Combo = paste(VariantType, PresentAbsentHigh, sep="_"))
AllVariants_SingleGene$FixedGeneNames = sapply(AllVariants_SingleGene$VariantLocus, function(x) str_replace(x,"\\'", "_"))
AllVariants_SingleGene$FixedGeneNames = sapply(AllVariants_SingleGene$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
AllVariants_SingleGene$FixedGeneNames = sapply(AllVariants_SingleGene$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
AllVariants_SingleGene$FixedGeneNames = sapply(AllVariants_SingleGene$FixedGeneNames, function(x) str_replace(x,":", "_"))
AllVariants_SingleGene$FixedGeneNames = sapply(AllVariants_SingleGene$FixedGeneNames, function(x) str_replace(x,"-", "_"))


AllVariants_SingleGene= AllVariants_SingleGene %>% left_join(GeneAnnotations,by="FixedGeneNames")
MostCommonByAnnotation = data.frame(table(AllVariants_SingleGene %>% select(Annotation, Comparison, Phenotype) %>% unique() %>% select(Annotation, Phenotype)))
MostCommonByAnnotation %>% arrange(-Freq)

# 7 comparisons in which staphylokinase associated with something annotated as a "Ig-like domain-containing protein" or 'transcriptional regulator'

AllVariants_SingleGene$TypePresenceCombo = paste(AllVariants_SingleGene$VariantType, AllVariants_SingleGene$PresentOrAbsentHigh, sep="_")

######################################
# Analyses of low-variant comparisons:
#####################################

# For each comparison, count how many total genes are affected

NumberGenesVaried = AllVariants_SingleGene %>% select(VariantLocus, Comparison) %>% unique() %>% group_by(Comparison) %>% summarize(NumGenesVaried = n())
AllVariants_SingleGene = AllVariants_SingleGene %>% left_join(NumberGenesVaried, by="Comparison")

View(AllVariants_HGTs)

AllVariants_HGTs %>% 
View(AllVariants_SingleGene %>% filter(NumGenesVaried <=10) %>% unique())

LowVariantComparisons = AllVariants_SingleGene %>% filter(NumGenesVaried <=10) %>% left_join(Present_In_Phage, by="FixedGeneNames") %>% left_join(Present_In_Plasmids, by="FixedGeneNames") %>% unique()
LowVariantComparisons[is.na(LowVariantComparisons)]<- 0 

LowVariantComparisons %>% arrange(Comparison, Comparison_Gene) %>% select(Gene,NumPhages, NumPlasmids )

View(LowVariantComparisons %>% select(Comparison_Gene, NumPhages, NumPlasmids) %>% unique())


OrderTypePresence = data.frame(ordertypes=1:6, TypePresenceCombo=c("AA_sub_STOPgained_absent","INDEL_FS_present",
                                                   "INDEL_FS_absent", "INDEL_inFrame_present",
                                                   "INDEL_inFrame_absent", "AA_sub_absent"))


View(WithPhagePlasmidInfo %>% select(Comparison_Gene, NumPlasmids, NumPhages))
WithPhagePlasmidInfo = (AllVariants_SingleGene %>% filter(NumVariants <= 10) %>% left_join(Present_In_Phage, by="FixedGeneNames") %>% left_join(Present_In_Plasmids, by="FixedGeneNames") %>% unique())
write.csv(WithPhagePlasmidInfo,"~/Documents/DataInputGithub/data/IntraPatient/Fewer10GenesVaried_PhagePlasmid.csv")



####################
# Biofilm barplot
####################
Biofilm_by_gene = AllVariants_SingleGene %>% filter(Phenotype=="Biofilm")

OccurrencesGeneBiofilm = table(Biofilm_by_gene %>% select(Comparison, VariantLocus) %>% unique() %>% select(VariantLocus)) 
OccurrencesNumPatientsBiofilm = table(Biofilm_by_gene %>% select(Patient, VariantLocus) %>% unique() %>% select(VariantLocus)) 
MoreThan1_genes = OccurrencesGeneBiofilm[OccurrencesGeneBiofilm>1]
MoreThan1_patient_Genes = OccurrencesNumPatientsBiofilm[OccurrencesNumPatientsBiofilm>1]
MoreThan1_patient_Annotation =table(Biofilm_by_gene %>% select(Annotation, Patient) %>% unique() %>% select(Annotation)) 
MoreThan1_patient_Annotation = MoreThan1_patient_Annotation[MoreThan1_patient_Annotation>1]

# no individual genes were consistent with biofilm in at least 2 patients
# however, 'hypothetical protein' and 'TIGR01741' were both annotations that came up twice:
OccurrencesAnnotationBiofilmMorethan1Patient = Biofilm_by_gene %>% filter(Annotation %in% names(MoreThan1_patient_Annotation)) %>% filter(Annotation!="hypothetical protein")

if(length(MoreThan1_patient_Genes)!=0){
  OccurrencesAnnotationBiofilm = table(Biofilm_by_gene %>% select(Comparison, Annotation) %>% unique() %>% select(Annotation)) 
  
  
  OccurrencesAnnotationBiofilm = table(Biofilm_by_gene %>% select(Comparison, Annotation) %>% unique() %>% select(Annotation)) 
  
  Genes_In_MoreThanOne_Comparison_Biofilm = names(sort(OccurrencesGeneBiofilm[OccurrencesGeneBiofilm>1]))
  
  MoreThan1GeneBiofilm = Biofilm_by_gene %>% filter(VariantLocus %in% Genes_In_MoreThanOne_Comparison_Biofilm)  
  AnnotationBarColorsBiofilm = colorScheme %>% filter(TypePresenceCombo %in% MoreThan1GeneBiofilm$TypePresenceCombo) 
  MoreThan1biofilmGenePlot = ggplot(MoreThan1GeneBiofilm, aes(x=VariantLocus, fill=TypePresenceCombo)) + geom_bar() + scale_fill_manual(values=AnnotationBarColorsBiofilm$HexVal) + coord_flip() + theme_classic()
  MoreThan1biofilmGenePlot$data$VariantLocus = factor(MoreThan1biofilmGenePlot$data$VariantLocus,levels=Genes_In_MoreThanOne_Comparison_Biofilm)
  
  
  BiofilmGeneSummary = MoreThan1GeneBiofilm %>% group_by(VariantLocus) %>%
    summarize(VariantLocus = VariantLocus,Annotation = Annotation, MeanNSvariants = mean(NumVariants),
              MeanNS_sequenceVariants = mean(NumSequenceVariants),NumComparisons=length(unique(Comparison)),
              NumPatients = length(unique(Patient)), CCs = paste(unique(CC), sep=";")) %>% 
    unique()
  write.csv(BiofilmGeneSummary, file="~/Documents/DataInputGithub/data/IntraPatient/Biofilm_MoreThanOneComparison_Genes.csv")
  
  
}



####################
# Siderophore barplot
####################
Siderophore_by_gene = AllVariants_SingleGene %>% filter(Phenotype=="Siderophore")

# different unique comparisons containing a variant for each gene
OccurrencesGeneSiderophore = table(Siderophore_by_gene %>% select(Comparison, VariantLocus) %>% unique() %>% select(VariantLocus))
bypatient = table(Siderophore_by_gene %>% select(Patient, VariantLocus) %>% unique() %>% select(VariantLocus))
moreThan1Patient = names(sort(bypatient[bypatient>1]))
Genes_In_MoreThanOne_Comparison_Siderophore = names(sort(OccurrencesGeneSiderophore[OccurrencesGeneSiderophore>1]))

# Plot the genes with highest # of comparisons found to have variants in that gene 
MoreThan1GeneSid = Siderophore_by_gene %>% filter(VariantLocus %in% moreThan1Patient)  

# Gene presence/absence isn't going to be duplicated for a given comparison
MoreThan1GeneSid_GeneLevel = MoreThan1GeneSid %>% filter(VariantType=="Gene")
MoreThan1GeneSid_SequenceLevel = MoreThan1GeneSid %>% filter(VariantType!="Gene")

MoreThan1GeneSid_Deduplicated = MoreThan1GeneSid_GeneLevel

for(gene in unique(MoreThan1GeneSid_SequenceLevel$Comparison_Gene)){
  subset_gene = MoreThan1GeneSid_SequenceLevel %>% filter(Comparison_Gene==gene)
  
  subset_gene = subset_gene %>% left_join(OrderTypePresence,by="TypePresenceCombo") %>% arrange(ordertypes) %>% select(-ordertypes)

  MoreThan1GeneSid_Deduplicated = rbind(MoreThan1GeneSid_Deduplicated, subset_gene[1,])
}

colorsToUse = MoreThan1GeneSid_Deduplicated %>% left_join(colorScheme, by="TypePresenceCombo") %>% arrange(TypePresenceCombo)
colorsToUse = unique(colorsToUse$Hex)


SiderophoreGeneSummary = MoreThan1GeneSid_Deduplicated %>% group_by(VariantLocus) %>%
  summarize(VariantLocus = VariantLocus,FixedGeneNames = FixedGeneNames,Annotation = Annotation, MeanGenesVaried = mean(NumGenesVaried),NumComparisons=length(unique(Comparison)),
            NumPatients = length(unique(Patient)), CCs = paste(unique(CC), collapse=";")) %>% arrange(NumComparisons, VariantLocus) %>% left_join(Present_In_Phage, by="FixedGeneNames") %>% left_join(Present_In_Plasmids, by="FixedGeneNames") %>% 
  unique()


MoreThan1siderophoreGene = ggplot(MoreThan1GeneSid_Deduplicated, aes(x=VariantLocus, fill=TypePresenceCombo)) + geom_bar() + scale_fill_manual(values=colorsToUse) + coord_flip() + theme_classic()
MoreThan1siderophoreGene$data$VariantLocus = factor(MoreThan1siderophoreGene$data$VariantLocus,levels=unique(SiderophoreGeneSummary$VariantLocus))

SiderophoreGeneSummary = data.frame(apply(SiderophoreGeneSummary, 2, rev))
SiderophoreGeneSummary$NumPhages = if_else(is.na(SiderophoreGeneSummary$NumPhages), 0, 1)
SiderophoreGeneSummary$NumPlasmids = if_else(is.na(SiderophoreGeneSummary$NumPlasmids), 0, 1)

write.csv(SiderophoreGeneSummary, file="~/Documents/DataInputGithub/data/IntraPatient/Siderophore_MoreThanOneComparison_Genes.csv")
ggsave(MoreThan1siderophoreGene, file="~/Documents/Saureus_Genomics_Paper/Figure4_SiderophoreConvergence.pdf", width=4.2,height=2)

SiderophoreGeneSummaryExpanded = MoreThan1GeneSid_Deduplicated %>% group_by(VariantLocus) %>%
  summarize(VariantLocus = VariantLocus,Annotation = Annotation, MeanGenesVaried = mean(NumGenesVaried),NumComparisons=length(unique(Comparison)),
            NumPatients = length(unique(Patient)), CCs = paste(unique(CC), collapse=";"),Patients = paste(unique(Patient), collapse=";")) %>% arrange(NumComparisons, VariantLocus) %>% 
  unique()


MoreThan1GeneSid_Deduplicated$patient_CC = paste(MoreThan1GeneSid_Deduplicated$Patient, MoreThan1GeneSid_Deduplicated$CC, sep="_")
SiderophoreGeneSummaryPlot = MoreThan1GeneSid_Deduplicated %>% group_by(VariantLocus) %>% select(patient_CC,Annotation,PresentOrAbsentHigh,VariantLocus, Comparison, CC) %>% unique()

SiderophoreGeneSummaryPlot = SiderophoreGeneSummaryPlot %>% arrange(VariantLocus, CC)


SiderophoreSummaryDF = data.frame()
for(gene in unique(SiderophoreGeneSummaryPlot$VariantLocus)){
  subDF=  SiderophoreGeneSummaryPlot %>% filter(VariantLocus==gene)
  subDF$Count = 1:nrow(subDF)
  SiderophoreSummaryDF = rbind(SiderophoreSummaryDF,subDF )
}
SiderophoreSummaryDF$CCLabel = SiderophoreSummaryDF$CC

siderophorecolors = unique((SiderophoreSummaryDF %>% left_join(cc_colors, by="CCLabel") %>% arrange(CC))$hexval)

siderophore_cc_count = ggplot(SiderophoreSummaryDF, aes(x=VariantLocus,y=factor(Count), fill=CC,color=CC)) +
  geom_star(starshape=13,size=3) + scale_fill_manual(values=siderophorecolors) + scale_color_manual(values=siderophorecolors) +
  theme_classic() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") + coord_flip()
siderophore_cc_count$data$VariantLocus = factor(siderophore_cc_count$data$VariantLocus, levels=rev(unique(SiderophoreGeneSummary$VariantLocus)))


ArrangedWithCount = cowplot::plot_grid(MoreThan1siderophoreGene, siderophore_cc_count, align="h",axis="l", rel_widths = c(3,2))
ggsave(ArrangedWithCount, file="~/Documents/Saureus_Genomics_Paper/Figure4_SiderophoreCCs.pdf",height=2, width=6)


SidLocus = sapply(SiderophoreSummaryDF$VariantLocus, function(x) str_replace(x, ":", "_"))
genepresence_phage %>% filter(FixedGeneNames %in% SidLocus) %>% select(FixedGeneNames,NumPhages)

genepresence_plasmids %>% filter(FixedGeneNames %in% SidLocus) %>% select(FixedGeneNames,NumPlasmids)

genepresence_plasmids %>% filter(FixedGeneNames=="group_237")
# patient 106
###############
MoreThan1GeneSid_Deduplicated %>% filter(Patient==106)
ComparatorGenomes106 = read.csv("~/Documents/DataInputGithub/data/IntraPatient/Distances/Comparisons_patient_106_CC8_Siderophore_1.csv")



########################
# Staphylokinase barplot
########################
Staphylokinase_by_gene = AllVariants_SingleGene %>% filter(Phenotype=="Staphylokinase")

bypatient = table(Staphylokinase_by_gene %>% select(Patient, VariantLocus) %>% unique() %>% select(VariantLocus))

bypatient = names(sort(bypatient[bypatient>2]))
OccurrencesGeneStaphylokinase = table(Staphylokinase_by_gene %>% select(Comparison, VariantLocus) %>% unique() %>% select(VariantLocus)) 
OccurrencesGeneStaphylokinase = OccurrencesGeneStaphylokinase[OccurrencesGeneStaphylokinase>2]


MoreThan1GeneStaphylokinase = Staphylokinase_by_gene %>% filter(VariantLocus %in% bypatient)

# choose representative of every gene's sequence variants 
MoreThan1GeneSak_GeneLevel = MoreThan1GeneStaphylokinase %>% filter(VariantType=="Gene")
MoreThan1GeneSak_SequenceLevel = MoreThan1GeneStaphylokinase %>% filter(VariantType!="Gene")
MoreThan1GeneSak_Deduplicated = MoreThan1GeneSak_GeneLevel
for(gene in unique(MoreThan1GeneSak_SequenceLevel$Comparison_Gene)){
  subset_gene = MoreThan1GeneSak_SequenceLevel %>% filter(Comparison_Gene==gene)
  
  subset_gene = subset_gene %>% left_join(OrderTypePresence,by="TypePresenceCombo") %>% arrange(ordertypes) %>% select(-ordertypes)
  
  MoreThan1GeneSak_Deduplicated = rbind(MoreThan1GeneSak_Deduplicated, subset_gene[1,])
}


AnnotationBarColorsSak = MoreThan1GeneSak_Deduplicated %>% left_join(colorScheme, by="TypePresenceCombo") %>% arrange(TypePresenceCombo)
colorsToUse = unique(AnnotationBarColorsSak$Hex)

KinaseGeneSummary = MoreThan1GeneStaphylokinase %>% group_by(VariantLocus) %>%
  summarize(VariantLocus = VariantLocus, FixedGeneNames = FixedGeneNames, Annotation = Annotation, MeanGenesVaried = mean(NumGenesVaried),NumComparisons=length(unique(Comparison)),
            NumPatients = length(unique(Patient)), CCs = paste(unique(CC), collapse=";")) %>% arrange(NumComparisons, VariantLocus) %>%left_join(Present_In_Phage, by="FixedGeneNames") %>% 
  left_join(Present_In_Plasmids, by="FixedGeneNames") %>% unique()


MoreThan1SakGenePlot = ggplot(MoreThan1GeneSak_Deduplicated, aes(x=VariantLocus, fill=TypePresenceCombo)) + geom_bar() + scale_fill_manual(values=colorsToUse) + coord_flip() + theme_classic()
MoreThan1SakGenePlot$data$VariantLocus = factor(MoreThan1SakGenePlot$data$VariantLocus,levels=(unique(KinaseGeneSummary$VariantLocus)))

KinaseGeneSummary =data.frame(apply((KinaseGeneSummary), 2, rev))

write.csv(KinaseGeneSummary, file="~/Documents/DataInputGithub/data/IntraPatient/Staphylokinase_MoreThanTwoPatient_Genes.csv")
ggsave(MoreThan1SakGenePlot, file="~/Documents/Saureus_Genomics_Paper/Figure4_StaphylokinaseConvergence.pdf", width=5)

unique(MoreThan1SakGenePlot$data$Patient)


KinaseGeneSummaryExpanded = MoreThan1GeneStaphylokinase %>% group_by(VariantLocus) %>%
  summarize(VariantLocus = VariantLocus, Annotation = Annotation, MeanGenesVaried = mean(NumGenesVaried),NumComparisons=length(unique(Comparison)),
            NumPatients = length(unique(Patient)), CCs = paste(unique(CC), collapse=";"), Patients = paste(unique(Patient), collapse=";")) %>% arrange(NumComparisons, VariantLocus) %>% unique()
KinaseGeneSummaryExpanded =data.frame(apply((KinaseGeneSummaryExpanded), 2, rev))




MoreThan1GeneSak_Deduplicated$patient_CC = paste(MoreThan1GeneSak_Deduplicated$Patient, MoreThan1GeneSak_Deduplicated$CC, sep="_")
StaphylokinaseGeneSummaryPlot = MoreThan1GeneSak_Deduplicated %>% group_by(VariantLocus) %>% select(patient_CC,Annotation,PresentOrAbsentHigh,VariantLocus, Comparison, CC) %>% unique()

StaphylokinaseGeneSummaryPlot = StaphylokinaseGeneSummaryPlot %>% arrange(VariantLocus, CC)


StaphylokinaseSummaryDF = data.frame()
for(gene in unique(StaphylokinaseGeneSummaryPlot$VariantLocus)){
  subDF=  StaphylokinaseGeneSummaryPlot %>% filter(VariantLocus==gene)
  subDF$Count = 1:nrow(subDF)
  StaphylokinaseSummaryDF = rbind(StaphylokinaseSummaryDF,subDF )
}
StaphylokinaseSummaryDF$CCLabel = StaphylokinaseSummaryDF$CC

kinasecolors = unique((StaphylokinaseSummaryDF %>% left_join(cc_colors, by="CCLabel") %>% arrange(CC))$hexval)

staphylokinase_cc_count = ggplot(StaphylokinaseSummaryDF, aes(x=VariantLocus,y=factor(Count), fill=CC,color=CC)) +
  geom_star(starshape=13,size=3) + scale_fill_manual(values=kinasecolors) + scale_color_manual(values=kinasecolors) +
  theme_classic() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") + coord_flip()
staphylokinase_cc_count$data$VariantLocus = factor(staphylokinase_cc_count$data$VariantLocus, levels=rev(unique(KinaseGeneSummary$VariantLocus)))


ArrangedWithCountSak = cowplot::plot_grid(MoreThan1SakGenePlot, staphylokinase_cc_count, align="h",axis="l", rel_widths = c(3,2))
ggsave(ArrangedWithCountSak, file="~/Documents/Saureus_Genomics_Paper/Figure4_StaphylokinaseCCs.pdf",height=5, width=6)



SakVarLocus=sapply(StaphylokinaseSummaryDF$VariantLocus, function(x) str_replace(x,":","_"))
genepresence_phage %>% filter(FixedGeneNames %in% SakVarLocus) %>% select(FixedGeneNames,NumPhages)
genepresence_plasmids %>% filter(FixedGeneNames %in% SakVarLocus) %>% select(FixedGeneNames,NumPlasmids)


########################
# Staphyloxanthin barplot
########################
Staphyloxanthin_by_gene = AllVariants_SingleGene %>% filter(Phenotype=="Staphyloxanthin")
OccurrencesGeneStaphyloxanthin= table(Staphyloxanthin_by_gene %>% select(Comparison, VariantLocus) %>% unique() %>% select(VariantLocus)) 

OccurrencesPatientStaphyloxanthin= table(Staphyloxanthin_by_gene %>% select(Patient, VariantLocus) %>% unique() %>% select(VariantLocus)) 
MoreThan1STX = OccurrencesGeneStaphyloxanthin[OccurrencesGeneStaphyloxanthin>1]

MoreThan1PatientStaphyloxanthin=names(OccurrencesPatientStaphyloxanthin[OccurrencesPatientStaphyloxanthin>1])

MoreThan1GeneSTX = Staphyloxanthin_by_gene %>% filter(VariantLocus %in% MoreThan1PatientStaphyloxanthin)

MoreThan1GeneSTX_GeneLevel = MoreThan1GeneSTX %>% filter(VariantType=="Gene")

# Only one sequence-level variant for STX
MoreThan1GeneSTX_SequenceLevel = MoreThan1GeneSTX %>% filter(VariantType!="Gene")
MoreThan1GeneSTX_Deduplicated = MoreThan1GeneSTX


AnnotationBarColorSTX = MoreThan1GeneSTX %>% left_join(colorScheme,by="TypePresenceCombo") %>% arrange(TypePresenceCombo) # filter(TypePresenceCombo %in% MoreThan1GeneSTX$TypePresenceCombo) 
colorsSTX = unique(AnnotationBarColorSTX$Hex)


STXGeneSummary = MoreThan1GeneSTX %>% group_by(VariantLocus) %>%
  summarize(FixedGeneNames=FixedGeneNames,VariantLocus = VariantLocus, Annotation = Annotation, MeanGenesVaried = mean(NumGenesVaried),NumComparisons=length(unique(Comparison)),
            NumPatients = length(unique(Patient)), CCs = paste(unique(CC), collapse=";")) %>% arrange(NumComparisons, VariantLocus) %>% left_join(Present_In_Phage, by="FixedGeneNames") %>% 
  left_join(Present_In_Plasmids, by="FixedGeneNames") %>% unique() 


MoreThan1STXGenePlot = ggplot(MoreThan1GeneSTX, aes(x=VariantLocus, fill=TypePresenceCombo)) + geom_bar() + scale_fill_manual(values=colorsSTX) + coord_flip() + theme_classic() + ylim(0,4)
MoreThan1STXGenePlot$data$VariantLocus = factor(MoreThan1STXGenePlot$data$VariantLocus,levels=(unique(STXGeneSummary$VariantLocus)))

STXGeneSummary = data.frame(apply(STXGeneSummary, 2, rev))
write.csv(STXGeneSummary, file="~/Documents/DataInputGithub/data/IntraPatient/Staphyloxanthin_MoreThanOneComparison_Genes.csv")

ggsave(MoreThan1STXGenePlot, file="~/Documents/Saureus_Genomics_Paper/Figure4_STXConvergence.pdf")
# Biofilm 201 clusters 1:3




MoreThan1GeneSTX_Deduplicated$patient_CC = paste(MoreThan1GeneSTX_Deduplicated$Patient, MoreThan1GeneSTX_Deduplicated$CC, sep="_")
STXGeneSummaryPlot = MoreThan1GeneSTX_Deduplicated %>% group_by(VariantLocus) %>% select(patient_CC,Annotation,PresentOrAbsentHigh,VariantLocus, Comparison, CC) %>% unique()

STXGeneSummaryPlot = STXGeneSummaryPlot %>% arrange(VariantLocus, CC)


STX_SummaryDF = data.frame()
for(gene in unique(STXGeneSummaryPlot$VariantLocus)){
  subDF=  STXGeneSummaryPlot %>% filter(VariantLocus==gene)
  subDF$Count = 1:nrow(subDF)
  STX_SummaryDF = rbind(STX_SummaryDF,subDF )
}
STX_SummaryDF$CCLabel = STX_SummaryDF$CC

STXcolors = unique((STX_SummaryDF %>% left_join(cc_colors, by="CCLabel") %>% arrange(CC))$hexval)

STX_cc_count = ggplot(STX_SummaryDF, aes(x=VariantLocus,y=factor(Count), fill=CC,color=CC)) +
  geom_star(starshape=13,size=3) + scale_fill_manual(values=STXcolors) + scale_color_manual(values=STXcolors) +
  theme_classic() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") + coord_flip()
STX_cc_count$data$VariantLocus = factor(STX_cc_count$data$VariantLocus, levels=rev(unique(STXGeneSummary$VariantLocus)))


ArrangedWithCounSTX = cowplot::plot_grid(MoreThan1STXGenePlot, STX_cc_count, align="h",axis="l", rel_widths = c(3,2))
ggsave(ArrangedWithCounSTX, file="~/Documents/Saureus_Genomics_Paper/Figure4_StxCCs.pdf",height=2.5, width=5.2)

STX_VariantLocus = sapply(STX_SummaryDF$VariantLocus, function(x) str_replace(x,":", "_"))
genepresence_phage %>% filter(FixedGeneNames %in% STX_VariantLocus) %>% select(FixedGeneNames,NumPhages)
genepresence_plasmids %>% filter(FixedGeneNames %in% STX_VariantLocus) %>% select(FixedGeneNames,NumPlasmids)


# Some examples to discuss
############################

# AC333
#
genomes162 = c("DORN1413", "DORN1399", "DORN1416", "DORN1410")

pat162_sak = AllVariants_SingleGene %>% filter(Comparison=="Staphylokinase_CC30_162_1_2")
In_AC33 = (genepresence_phage_plasmid %>% filter(AC333==1))
In_AC33$FixedGeneNames = sapply(In_AC33$X, function(x) str_replace(x,"\\'", "_"))
In_AC33$FixedGeneNames = sapply(In_AC33$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
In_AC33$FixedGeneNames = sapply(In_AC33$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
In_AC33$FixedGeneNames = sapply(In_AC33$FixedGeneNames, function(x) str_replace(x,":", "_"))
In_AC33$FixedGeneNames = sapply(In_AC33$FixedGeneNames, function(x) str_replace(x,"-", "_"))

setdiff(pat162_sak$FixedGeneNames, In_AC33$FixedGeneNames )
GeneAnnotations %>% filter(Gene %in%In_AC33$X )
phagepresence_isolates %>% filter(X %in% genomes162)

# all contain Phage61
# 


# Sak in Patient 141
pat141_sak = AllVariants_SingleGene %>% filter(Comparison=="Staphylokinase_CC1_141_1_3")
pat141_sak2_3 = AllVariants_SingleGene %>% filter(Comparison=="Staphylokinase_CC1_141_2_3")

In_52 = (genepresence_phage_plasmid %>% filter(Phage52==1))
In_52$FixedGeneNames = sapply(In_52$X, function(x) str_replace(x,"\\'", "_"))
In_52$FixedGeneNames = sapply(In_52$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
In_52$FixedGeneNames = sapply(In_52$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
In_52$FixedGeneNames = sapply(In_52$FixedGeneNames, function(x) str_replace(x,":", "_"))
In_52$FixedGeneNames = sapply(In_52$FixedGeneNames, function(x) str_replace(x,"-", "_"))



setdiff(pat141_sak$FixedGeneNames, In_52$FixedGeneNames)

setdiff(pat141_sak2_3$FixedGeneNames, In_52$FixedGeneNames)


GeneAnnotations %>% filter(Gene %in%In_AC33$X )
phagepresence_isolates %>% filter(X %in% genomes162)

GeneAnnotations162 %>% filter(GeneAnnotations=="sak")
