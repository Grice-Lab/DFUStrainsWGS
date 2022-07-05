
# Amy Campbell
# 04/2022
# Differential expression analysis for DORN925 and DORN1088 exposed to H202 in vitro 
library(dplyr)
library(Rsubread)
library(Rsamtools)
library(stringr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(tidyr)
library(viridis)

samplemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/MappingRNASeq.txt", sep="\t", header=F)
colnames(samplemap) = c("Sample", "SampleTreatment", "Condition")
genemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/GeneNamesPreliminary.csv", sep=',')
load("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNACounts.rda")

# make GeneAdjust column where if no gene name available, list RefSeq ID with homology as ID'd by PGAP 
########################################################################################################
genemap = genemap %>% mutate(GeneAdjust= if_else( (is.na(Gene) | Gene==""), RefSeqID, Gene))

# Manipulate sample ID mapping file and make plot of Sample vs. Total Reads assigned
#####################################################################################
CountsBySample=data.frame((colSums(MyCounts$counts)))
CountsBySample$Sample=sapply(row.names(CountsBySample), function(x) str_split(x, pattern="_")[[1]][1])
colnames(CountsBySample) = c("TotalAssignedReads", "Sample")

CountsBySample = samplemap %>% left_join(CountsBySample, by="Sample")
CountsBySample[is.na(CountsBySample)] <- 0
CountsBySample$Strain = sapply(CountsBySample$SampleTreatment, function(x) if_else(grepl(x, pattern="ARM72"),"DORN925(high)","DORN1088(low)") )

CountsBySample$ConditionStrain = paste(CountsBySample$Strain, CountsBySample$Condition, sep=" : ")
CountsBySample$ConditionStrain = factor(CountsBySample$ConditionStrain )
TotalCountsPlot = ggplot(CountsBySample, aes(x=Sample, y=TotalAssignedReads, fill=ConditionStrain)) + geom_bar(stat="identity", color="black")
TotalCountsPlot$data$Sample = factor(TotalCountsPlot$data$Sample , levels=(TotalCountsPlot$data %>% arrange(ConditionStrain))$Sample)

TotalCountsPlot$data$Sample = factor(TotalCountsPlot$data$Sample , levels=(TotalCountsPlot$data %>% arrange(Strain, Condition))$Sample)

# Low_control, Low_Treatment, High_control, High_Treatment
treatmentcombos = c("#BCD2E8", "#2E5984", "#FCF787", "#EAAA00")
TotalCountsPlot=TotalCountsPlot + scale_fill_manual(values=treatmentcombos) + scale_y_continuous(limits=c(0, 25000000), breaks=seq(0, 25000000, 1000000))
ggsave(TotalCountsPlot,file= "/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/ReadCountsAssigned.png")

SampleMapping = CountsBySample %>% filter(Sample!="ARM6")
CountsBySample$TotalAssignedReads = NULL

# Read counts data into DEseq
##############################
columndataDE = SampleMapping %>% select(Condition, Strain)
columndataDE$Condition= factor(columndataDE$Condition)

columndataDE$Strain = sapply(columndataDE$Strain, function(x) str_split(x, pattern='\\(')[[1]][1])
columndataDE$Strain= factor(columndataDE$Strain)

rownames(columndataDE) = SampleMapping$Sample

CountsMatrix = MyCounts$counts
colnames(CountsMatrix) = sapply(colnames(CountsMatrix), function(x) str_split(x, pattern="_")[[1]][1])
CountsMatrix = CountsMatrix[,row.names(columndataDE)]

# "group" is the condition_strain (so ctrl_DORN925, treatment_925, etc.)
columndataDE$group = factor(paste(columndataDE$Condition, columndataDE$Strain, sep="_"))

DEData =DESeqDataSetFromMatrix(countData=CountsMatrix, 
                               colData=columndataDE,
                               design= ~group)
myDESeqObj = DESeq(DEData)


# PCA on VST-transformed
########################
myDESeqObjVST = vst(myDESeqObj)
plotPCA(myDESeqObjVST,intgroup=c("group")) + scale_color_manual()

# Sanity check:sak should be 0 counts in DORN1088 samples
##########################################################
CountsMatrix[c("cds-pgaptmp_000025"),]
# 0 in all of what are supposed to be 1088 genomes so that's good. (ARMs 3,4, 7,8, 11, 12)


###############################################
###############################################
# (1) How does DORN925 respond to H202 stress?
###############################################
###############################################

# This should show results for 
# (expression after treatment in DORN925)/(expression before treatment in DORN925) AKA log2(treatment_ex/control_ex)
results_DORN925_Stress <- results(myDESeqObj, contrast=c("group", "treatment_DORN925","ctrl_DORN925"))
results_DORN925_StressDF = (data.frame(results_DORN925_Stress))

# Add Uniprot, gene names, GOs
results_DORN925_StressDF$AnnotID = row.names(results_DORN925_StressDF)
results_DORN925_StressDF = results_DORN925_StressDF %>% left_join(genemap, by="AnnotID")


# Make some variables for plotting & filtering
##############################################

# DE: whether log-fold change > ±.6 (~ 1.5X fold change) and whether FDR- adjusted p value < .05
results_DORN925_StressDF$DE = if_else( abs(results_DORN925_StressDF$log2FoldChange) < .6  | results_DORN925_StressDF$padj >= .05 | is.na(results_DORN925_StressDF$padj), FALSE, TRUE)

# DirectionDE: whether DE is true and whether up-regulated or down-regulated if so
results_DORN925_StressDF = results_DORN925_StressDF %>% mutate(DirectionDE = case_when( log2FoldChange > 0 & DE==TRUE ~ "Up",
                                                                                        log2FoldChange <0 & DE== TRUE ~ "Down", 
                                                                                        TRUE ~ "None"))

# Label: "Up" / "Down" if the total log2foldchange >2 (4x-fold change) and DE is true; NA otherwise (so as not to crowd the volcano plot)
results_DORN925_StressDF = results_DORN925_StressDF %>% mutate(Label = if_else( ((DirectionDE == "Up" | DirectionDE == "Down" ) & abs(log2FoldChange) >2 ), GeneAdjust,""))
results_DORN925_StressDF$Label[results_DORN925_StressDF$Label==""] <- NA



resultsStress925 = results_DORN925_StressDF %>% select(log2FoldChange, padj, AnnotID, RefSeqID, UniprotID, Gene)
write.csv(resultsStress925, file="data/RNASeqDataSets/TotalGenesStressVsControl925.csv")

# Volcano plot for 925 -- everything that fits with a >4X fold change and < .05 p-value labeled 
################################################################################################
VolcPlot925 = ggplot(data=results_DORN925_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=Label)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel() + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN925\nin Response to H202 Stress") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5))

# MA plot for 925
####################
labels=results_DORN925_StressDF$Label[!is.na(results_DORN925_StressDF$Label)]
MAplot925 = ggmaplot(results_DORN925_StressDF,size=1, genenames = results_DORN925_StressDF$GeneAdjust, label.select=labels) + ggtitle("Differentially expressed genes in DORN925\nin Response to H202 Stress")
MAplot925= MAplot925 + theme(plot.title=element_text(face="bold", size=20, hjust=.5))

ggsave(MAplot925, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotDORN925Stress.png", height=10, width=20)


#############################################################
#############################################################
# (2) How does DORN1088 (low xanthin) respond to H202 stress? 
#############################################################
#############################################################


## (expression after treatment in DORN1088)/(expression before treatment in DORN1088) AKA log2(treatment_ex/control_ex)
results_DORN1088_Stress <- results(myDESeqObj, contrast=c("group", "treatment_DORN1088","ctrl_DORN1088"))
results_DORN1088_StressDF = (data.frame(results_DORN1088_Stress))

# Add Uniprot, gene names, GOs
results_DORN1088_StressDF$AnnotID = row.names(results_DORN1088_StressDF)
results_DORN1088_StressDF = results_DORN1088_StressDF %>% left_join(genemap, by="AnnotID")

# Make some variables for plotting & filtering
############################################## 

# DE: whether log-fold change > ±.6 (~ 1.5X fold change) and whether FDR- adjusted p value < .05
results_DORN1088_StressDF$DE = if_else( abs(results_DORN1088_StressDF$log2FoldChange) < .6  | results_DORN1088_StressDF$padj >= .05 | is.na(results_DORN1088_StressDF$padj), FALSE, TRUE)

# DirectionDE: whether DE is true and whether up-regulated or down-regulated if so
results_DORN1088_StressDF = results_DORN1088_StressDF %>% mutate(DirectionDE = case_when( log2FoldChange > 0 & DE==TRUE ~ "Up",
                                                                                          log2FoldChange <0 & DE== TRUE ~ "Down", 
                                                                                          TRUE ~ "None"))
# only label those with 4X and higher (log2Foldchange>=2)
results_DORN1088_StressDF = results_DORN1088_StressDF %>% mutate(Label = if_else( ((DirectionDE == "Up" | DirectionDE == "Down" ) & abs(log2FoldChange) >2 ), GeneAdjust, ""))
results_DORN1088_StressDF$Label[results_DORN1088_StressDF$Label==""] <- NA


resultsStress1088 = results_DORN1088_StressDF %>% select(log2FoldChange, padj, AnnotID, RefSeqID, UniprotID, Gene)
write.csv(resultsStress1088, file="data/RNASeqDataSets/TotalGenesStressVsControl1088.csv")

VolcPlot1088 = ggplot(data=results_DORN1088_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=Label)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel() + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN1088\nin Response to H202 Stress") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5)) +



# MA plot1088
#######################
mylabels = results_DORN1088_StressDF$Label[!is.na(results_DORN1088_StressDF$Label)]
MAplot1088 = ggmaplot(results_DORN1088_StressDF,size=1, genenames = results_DORN1088_StressDF$Label, label.select=mylabels) + ggtitle("Differentially expressed genes in DORN1088\nin Response to H202 Stress")

MAplot1088= MAplot1088 + theme(plot.title=element_text(face="bold", size=20, hjust=.5))

ggsave(MAplot1088, file="data/RNASeqPlots/MAPlotDORN1088Stress.png", height=10, width=20)


gridarrange= gridExtra::grid.arrange(VolcPlot1088, VolcPlot925, ncol=2)
ggsave(gridarrange, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/SideBySide.pdf", width=30, height=15)

# Upregulated genes in DORN1088 in response to H202 stress
UP_1088 = results_DORN1088_StressDF %>% filter(DirectionDE=="Up")
write.csv(UP_1088 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/UpDE_ARM81.csv")

# Downregulated genes in DORN1088 in response to H202 stress 
DOWN_1088 = results_DORN1088_StressDF %>% filter(DirectionDE=="Down")
write.csv(DOWN_1088 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DownDE_ARM81.csv")

# Upregulated genes in DORN925 in response to H202 stress
UP_925 = results_DORN925_StressDF %>% filter(DirectionDE=="Up")
write.csv(UP_925 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/UpDE_ARM72.csv")

# Downregulated genes in DORN925 in response to H202 stress 
DOWN_925 = results_DORN925_StressDF %>% filter(DirectionDE=="Down")
write.csv(DOWN_925 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DownDE_ARM72.csv")


###################################################################################
###################################################################################
# (3)Which genes are uniquely Up/Down-regulated in response to H202 in each strain? 
###################################################################################
###################################################################################

UniqueUp925 = results_DORN925_StressDF %>% filter(AnnotID %in% setdiff(UP_925$AnnotID, UP_1088$AnnotID))
UniqueDown925 = results_DORN925_StressDF %>% filter( (AnnotID %in% setdiff(DOWN_925$AnnotID, DOWN_1088$AnnotID)) )
UniqueUp1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% setdiff(UP_1088$AnnotID,UP_925$AnnotID) )
UniqueDown1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% setdiff(DOWN_1088$AnnotID,DOWN_925$AnnotID) )

write.csv(UniqueUp925 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/UpDE_ARM72_Unique.csv")
write.csv(UniqueDown925 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DownDE_ARM72_Unique.csv")
write.csv(UniqueDown1088 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DownDE_ARM81_Unique.csv")
write.csv(UniqueUp1088 %>% select(log2FoldChange, pvalue, padj, AnnotID, RefSeqID, UniprotID, GeneAdjust, GO_terms, DirectionDE), file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/UpDE_ARM81_Unique.csv")


results_DORN925_StressDF = results_DORN925_StressDF %>% mutate(LabelUnique = if_else(DE==TRUE & (AnnotID %in%   append(UniqueUp925$AnnotID, UniqueDown925$AnnotID)),GeneAdjust, ""))
results_DORN925_StressDF$LabelUnique[results_DORN925_StressDF$LabelUnique ==""] <- NA

# Volcano plot but label only the genes that fit and are uniquely up or down-regulated in DORN925
VolcPlot925UniqueGr2 = ggplot(data=results_DORN925_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=LabelUnique)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel(color="black") + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN925\nin Response to H202 Stress (Unique)") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5))


results_DORN1088_StressDF = results_DORN1088_StressDF %>% mutate(LabelUnique = if_else(DE==TRUE & (AnnotID %in%   append(UniqueUp1088$AnnotID, UniqueDown1088$AnnotID)),GeneAdjust, ""))
results_DORN1088_StressDF$LabelUnique[results_DORN1088_StressDF$LabelUnique ==""] <- NA

# Volcano plot but label only the genes that fit and are uniquely up or down-regulated in DORN925
VolcPlot1088UniqueGr2 = ggplot(data=results_DORN1088_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=LabelUnique)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel(color="black") + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN1088\nin Response to H202 Stress(Unique)") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5))


arranged = gridExtra::grid.arrange(VolcPlot1088UniqueGr2,VolcPlot925UniqueGr2,ncol=2)
ggsave(arranged, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/SideBySideUnique.pdf", width=30, height=15)


###################################################################################
###################################################################################
# (4) Which genes are DE between the two strains under control conditions? H202?
###################################################################################
###################################################################################


results_genotype <- results(myDESeqObj, contrast=c("group", "ctrl_DORN925","ctrl_DORN1088"))

View(data.frame(results_genotype))

MapGenes = data.frame(AnnotID=row.names(results_genotype))
MapGenes = MapGenes %>% left_join(genemap, by="AnnotID") %>% unique()

MAplotStrain = ggmaplot(results_genotype,size=1, genenames = MapGenes$GeneAdjust)
MAplotStrain = MAplotStrain + ggtitle("Differentially expressed genes between \nDORN925:DORN1088 under control conditions") + theme(plot.title=element_text(face="bold", size=20, hjust=.5))
ggsave(MAplotStrain, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotControlStrains.png")

results_genotypeDF = data.frame(results_genotype) %>% select(padj,pvalue, log2FoldChange)
results_genotypeDF$AnnotID = row.names(results_genotypeDF)
results_genotypeDF = results_genotypeDF %>% left_join(genemap, by="AnnotID")
results_genotypeDF$DE = if_else( abs(results_genotypeDF$log2FoldChange) < .6  | results_genotypeDF$padj >= .05 | is.na(results_genotypeDF$padj), FALSE, TRUE)
write.csv(results_genotypeDF, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DORN925_to_DORN1088_Control.csv")


results_genotypeH202 <- results(myDESeqObj, contrast=c("group", "treatment_DORN925","treatment_DORN1088"))
results_genotypeH202DF = data.frame(results_genotypeH202) %>% select(padj,pvalue, log2FoldChange)
results_genotypeH202DF$DE = if_else( abs(results_genotypeH202DF$log2FoldChange) < .6  | results_genotypeH202DF$padj >= .05 | is.na(results_genotypeH202DF$padj), FALSE, TRUE)


MapGenesH202 = data.frame(AnnotID = row.names(results_genotypeH202DF))
MapGenesH202 = MapGenesH202 %>% left_join(genemap, by="AnnotID") %>% unique()

MAplotStrainH202 = ggmaplot(results_genotypeH202,size=1, genenames = MapGenesH202$GeneAdjust) +
  ggtitle("Differentially expressed genes between \nDORN925:DORN1088 after H202 exposure") + theme(plot.title=element_text(face="bold", size=20, hjust=.5))
ggsave(MAplotStrainH202, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotH202Strains.png")

results_genotypeH202DF$AnnotID = row.names(results_genotypeH202DF)
results_genotypeH202DF = results_genotypeH202DF %>% left_join(genemap, by="AnnotID")
write.csv(results_genotypeH202DF,"/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DORN925_to_DORN1088_H202.csv")


UP_Genotype_Control = results_genotypeDF %>% filter(DE==TRUE & log2FoldChange>0)
DOWN_Genotype_Control = results_genotypeDF %>% filter(DE==TRUE & log2FoldChange<0)

UP_Genotype_H202 = results_genotypeH202DF %>% filter(DE==TRUE & log2FoldChange>0)
DOWN_Genotype_H202 = results_genotypeH202DF %>% filter(DE==TRUE & log2FoldChange<0)

UniqueDOWNControl = DOWN_Genotype_Control %>% filter(AnnotID %in% setdiff(DOWN_Genotype_Control$AnnotID, DOWN_Genotype_H202$AnnotID))
UniqueDOWNH202 = DOWN_Genotype_H202 %>% filter(AnnotID %in% setdiff(DOWN_Genotype_H202$AnnotID, DOWN_Genotype_Control$AnnotID))
write.csv(UniqueDOWNControl,"/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DOWN_DORN925_to_DORN1088_UniqueControl.csv")
write.csv(UniqueDOWNH202,"/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/DOWN_DORN925_to_DORN1088_UniqueH202.csv")



UniqueUpControl = UP_Genotype_Control %>% filter(AnnotID %in% setdiff(UP_Genotype_Control$AnnotID, UP_Genotype_H202$AnnotID))
UniqueUpH202 = UP_Genotype_H202 %>% filter(AnnotID %in% setdiff(UP_Genotype_H202$AnnotID, UP_Genotype_Control$AnnotID))

write.csv(UniqueUpControl,"/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/UP_DORN925_to_DORN1088_UniqueControl.csv")


write.csv(UniqueUpH202,"/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqDataSets/UP_DORN925_to_DORN1088_UniqueH202.csv")




###########################################################################################################
##########################################################################################################
# (5) How are common stress genes DE inresponse to H202 in each strain? Between the strains at baseline?
########################################################################################################
##########################################################################################################


# Looking just at  genes of interest in each strain in response to stress
###########################################################################
stressGenes = read.csv("data/StressResponseGenesInclude.csv")
stressGenes = stressGenes %>% filter(! (GeneName %in% c("psmA1",
                                     "psmA2",
                                     "psmA3",
                                     "psmA4","hld")))
# Stress genes in DORN925 re: stress
results_DORN925_StressDFjoined = results_DORN925_StressDF %>% left_join(stressGenes, by="AnnotID")
StressResponsive925stress = results_DORN925_StressDFjoined %>% filter(!is.na(GeneName))
StressResponsive925stress = StressResponsive925stress %>% group_by(Grouping, OperonRegulonSubgrouping) %>% arrange(Grouping,OperonRegulonSubgrouping,  GeneName) %>% ungroup()
StressResponsive925stress = StressResponsive925stress %>% mutate(GeneAdjust= if_else(is.na(GeneName), GeneAdjust, GeneName))
StressResponsive925stress$comparison = "SA925_h2o2"
StressResponsive925 = StressResponsive925stress %>% select(GeneAdjust, comparison, log2FoldChange, padj)


order = StressResponsive925stress$GeneAdjust
# Stress genes in DORN1088 re: stress
results_DORN1088_Stressjoined = results_DORN1088_StressDF %>% left_join(stressGenes, by="AnnotID")
results_DORN1088stress = results_DORN1088_Stressjoined %>% filter(!is.na(GeneName))
results_DORN1088stress = results_DORN1088stress %>% group_by(Grouping, OperonRegulonSubgrouping) %>% arrange(Grouping,OperonRegulonSubgrouping,  GeneName) %>% ungroup()
StressResponsive1088stress = results_DORN1088stress %>% mutate(GeneAdjust= if_else(is.na(GeneName), GeneAdjust, GeneName))
StressResponsive1088stress$comparison = "SA1088_h2o2"
StressResponsive1088stress = StressResponsive1088stress %>% select(GeneAdjust, comparison, log2FoldChange, padj)

# Stress genes in 925 vs. 1088 under control conditions
results_Genotype_joined = results_genotypeDF %>% left_join(stressGenes, by="AnnotID")
results_Genotype_joined = results_Genotype_joined %>% filter(!is.na(GeneName))
results_Genotype_joined = results_Genotype_joined %>% group_by(Grouping, OperonRegulonSubgrouping) %>% arrange(Grouping,OperonRegulonSubgrouping,  GeneName) %>% ungroup()
StressResponsiveGenotype = results_Genotype_joined %>% mutate(GeneAdjust= if_else(is.na(GeneName), GeneAdjust, GeneName))
StressResponsiveGenotype$comparison = "SA925_1088_Control"
StressResponsiveGenotype = StressResponsiveGenotype %>% select(GeneAdjust, comparison, log2FoldChange, padj)


View(results_Genotype_joined)

StressGenesAll = rbind(StressResponsive925, StressResponsive1088stress)
StressGenesAll = rbind(StressGenesAll, StressResponsiveGenotype)

StressGenesMelted = StressGenesAll %>% reshape2::melt(id.vars=c("comparison", "GeneAdjust"))
StressGenesMeltedPvalues =  StressGenesMelted %>% filter(variable=="padj")
  
StressGenesMeltedLFCs =  StressGenesMelted %>% filter(variable=="log2FoldChange")

StressGenesMeltedPvalues = StressGenesMeltedPvalues %>%  mutate(sigLabel=case_when(value < .05 & value >= .01 ~"*",
                                                               value <.01 & value >= .001 ~"**",
                                                               value < .001 & value >= .0001 ~ "***",
                                                               value < .0001~ "****",

                                                               TRUE~ ""))
StressGenesMeltedPvalues
StressGenesMeltedLFCs = StressGenesMeltedLFCs %>% left_join(StressGenesMeltedPvalues %>% select(comparison, GeneAdjust, sigLabel), by=c("comparison", "GeneAdjust"))



#View(StressGenesMeltedLFCs)
testplot = ggplot(StressGenesMeltedLFCs, aes(x=comparison, y=GeneAdjust, fill=value)) + geom_tile() + scale_fill_viridis(option="plasma") + geom_text(aes(label=sigLabel), vjust=.5, color="white") + theme_classic() + labs(fill="log2FC") + ylab("Gene")
testplot$data$GeneAdjust = factor(testplot$data$GeneAdjust, levels=rev(order))
testplot$data$comparison = factor(testplot$data$comparison, levels=c("SA925_1088_Control", "SA925_h2o2", "SA1088_h2o2"))
testplot
testplot + scale_fill_viridis(option="inferno")

ggplot(StressGenesMeltedLFCs, aes(x=comparison, y=GeneAdjust, fill=value)) + geom_tile() + scale_fill_gradient(high="#FFFF00" , low="#800080") + geom_text(aes(label=sigLabel), vjust=.5, color="white") + theme_classic() + labs(fill="log2FC") + ylab("Gene")



#FFD014 60339B

#scale_fill_gradient2(high="red", mid="white", low="blue", midpoint = 0)

crtGenes = c("cds-pgaptmp_002177", "cds-pgaptmp_002178", "cds-pgaptmp_002179", "cds-pgaptmp_002180", "cds-pgaptmp_002181")
sigBGenes=c("cds-pgaptmp_002685", "cds-pgaptmp_002686", "cds-pgaptmp_002687", "cds-pgaptmp_002688")

perGenes = c("cds-pgaptmp_000116","cds-pgaptmp_000082","cds-pgaptmp_000114","cds-pgaptmp_000469", "cds-pgaptmp_000756","cds-pgaptmp_000957", "cds-pgaptmp_001682", "cds-pgaptmp_001683", "cds-pgaptmp_001275")
msrGenes = c("cds-pgaptmp_000661", "cds-pgaptmp_000662" )
clpGenes = c("cds-pgaptmp_001530", "cds-pgaptmp_001271", "cds-pgaptmp_001133", "cds-pgaptmp_000349", "cds-pgaptmp_002195", "cds-pgaptmp_001531", "cds-pgaptmp_001532")
sodGenes = c("cds-pgaptmp_000470", "cds-pgaptmp_001935")

phageGenesAll = read.csv("data/AnnotationsPhage925_PGAP_Patric.csv")
phageGenesAllList = sapply(phageGenesAll$PGAPid,function(x) (str_split(x, pattern="ID="))[[1]][2])
phageGenesAll$AnnotID = sapply(phageGenesAll$PGAPid,function(x) (str_split(x, pattern="ID="))[[1]][2] )

# SAOUHSC_01999 is bcp 
# "cds-pgaptmp_002613" = mrgA based on tblastn of b subtilis's mgrA 
# "cds-pgaptmp_001532" = mcsA based on tblastn of s auerus mcsA
results_DORN925_StressDF$GeneAdjust = if_else(results_DORN925_StressDF$GeneAdjust=="SAOUHSC_01999", "bcp", results_DORN925_StressDF$GeneAdjust)
results_DORN1088_StressDF$GeneAdjust = if_else(results_DORN1088_StressDF$GeneAdjust=="SAOUHSC_01999", "bcp", results_DORN1088_StressDF$GeneAdjust)
results_genotypeDF$GeneAdjust = if_else(results_genotypeDF$GeneAdjust=="SAOUHSC_01999", "bcp", results_genotypeDF$GeneAdjust)

results_DORN925_StressDF$GeneAdjust = if_else(results_DORN925_StressDF$AnnotID=="cds-pgaptmp_002613", "mgrA", results_DORN925_StressDF$GeneAdjust)
results_DORN1088_StressDF$GeneAdjust = if_else(results_DORN1088_StressDF$AnnotID=="cds-pgaptmp_002613", "mgrA", results_DORN1088_StressDF$GeneAdjust)
results_genotypeDF$GeneAdjust = if_else(results_genotypeDF$AnnotID=="cds-pgaptmp_002613", "mgrA", results_genotypeDF$GeneAdjust)



results_DORN925_StressDF$GeneAdjust = if_else(results_DORN925_StressDF$AnnotID=="cds-pgaptmp_001532", "mcsA", results_DORN925_StressDF$GeneAdjust)
results_DORN1088_StressDF$GeneAdjust = if_else(results_DORN1088_StressDF$AnnotID=="cds-pgaptmp_001532", "mcsA", results_DORN1088_StressDF$GeneAdjust)
results_genotypeDF$GeneAdjust = if_else(results_genotypeDF$AnnotID=="cds-pgaptmp_001532", "mcsA", results_genotypeDF$GeneAdjust)



crts925 = results_DORN925_StressDF %>% filter(AnnotID %in% crtGenes) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID)
crts1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% crtGenes) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID)
crtsGenotype = results_genotypeDF %>% filter(AnnotID %in% crtGenes) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID)



sig925 = results_DORN925_StressDF %>% filter(AnnotID %in% sigBGenes) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 
sig1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% sigBGenes) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 

write.csv(crts925, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/crt925.csv", row.names=F)
write.csv(sig925, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/sig925.csv", row.names=F)

write.csv(crts1088, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/crt1088.csv", row.names=F)
write.csv(sig1088, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/sig1088.csv", row.names=F)

sods925 = results_DORN925_StressDF %>% filter(AnnotID %in% c("cds-pgaptmp_000470", "cds-pgaptmp_001935")) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 
sods1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% c("cds-pgaptmp_000470", "cds-pgaptmp_001935")) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 

write.csv(sods925, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/sod925.csv", row.names=F)
write.csv(sods1088, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/sod1088.csv", row.names=F)

kat925 = results_DORN925_StressDF %>% filter(AnnotID %in% c("cds-pgaptmp_000756")) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 
kat1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% c("cds-pgaptmp_000756")) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 

per925 = results_DORN925_StressDF %>% filter(AnnotID %in% perGenes) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID)
per1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% perGenes) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID)

write.csv(per925, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/per925.csv", row.names=F)
write.csv(per1088, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/per1088.csv", row.names=F)

agr925 = results_DORN925_StressDF %>% filter(GeneAdjust %in% c("agrA", "agrB", "agrC", "agrD")) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 
agr1088 = results_DORN1088_StressDF %>% filter(GeneAdjust %in% c("agrA", "agrB", "agrC", "agrD")) %>% select(GeneAdjust, log2FoldChange, pvalue, padj, AnnotID) 
write.csv(agr925, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/agr925.csv", row.names=F)
write.csv(agr1088, quote=F, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/agr1088.csv", row.names=F)

phage925 = results_DORN925_StressDF %>% filter(AnnotID %in% phageGenesAllList)
View(phage925 %>% filter(DE=="TRUE") %>% left_join(phageGenesAll, by="AnnotID"))

Sig925 = rbind(crts925, sig925, agr925)
Sig1088 = rbind(crts1088, sig1088, agr1088)

Sig925 = rbind(crts925, sig925, agr925)
Sig1088 = rbind(crts1088, sig1088, agr1088)
Sig925$Strain="DORN925"
Sig1088$Strain="DORN1088"





StressGenes = c(perGenes,msrGenes, clpGenes, sodGenes )

Stress925 = results_DORN925_StressDF %>% filter(AnnotID %in% StressGenes)
Stress1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% StressGenes)
Stress925$Strain="DORN925"
Stress1088$Strain="DORN1088"
StressGenotype = results_genotypeDF %>% filter(AnnotID %in% StressGenes)


MeltedSigmaB = rbind(Sig1088, Sig925) %>% select(GeneAdjust, Strain, log2FoldChange, padj) %>% reshape2::melt(id.vars=c("Strain", "GeneAdjust"))
MeltedSigmaBPs = MeltedSigmaB %>% filter(variable=="padj")
MeltedSigmaB_FCs = MeltedSigmaB %>% filter(variable=="log2FoldChange")
MeltedSigmaBPs = MeltedSigmaBPs %>%  mutate(sigLabel=case_when(value < .05 & value >= .01 ~"*", 
                                                               value <.01 & value >= .001 ~"**",
                                                               value < .001 & value >= .0001 ~ "***", 
                                                               value < .0001~ "****",
                                                               
                                                               TRUE~ ""))
MeltedSigmaB_FCs = MeltedSigmaB_FCs %>% left_join(MeltedSigmaBPs %>% select(Strain, GeneAdjust, sigLabel), by=c("Strain", "GeneAdjust"))
ggplot(MeltedSigmaB_FCs, aes(x=Strain, y=GeneAdjust, fill=value)) + geom_tile() + scale_fill_gradient2(high="red", mid="white", low="blue", midpoint = 0)+ geom_text(aes(label=sigLabel)) + theme_classic() + labs(fill="log2FC after H202") + ylab("Gene")



MeltedStress = rbind(Stress925, Stress1088) %>% select(GeneAdjust, Strain, log2FoldChange, padj) %>% reshape2::melt(id.vars=c("Strain", "GeneAdjust")) 
MeltedStress$GeneAdjust = if_else(MeltedStress$GeneAdjust=="msrA2", "msrA", MeltedStress$GeneAdjust)

MeltedStressPs = MeltedStress %>% filter(variable=="padj")
MeltedStressPs = MeltedStressPs %>% mutate(sigLabel=case_when(value < .05 & value >= .01 ~"*", 
                                                              value <.01 & value >= .001 ~"**",
                                                              value < .001 & value >= .0001 ~ "***", 
                                                              value < .0001~ "****",
                                                              
                                                              TRUE~ ""))

MeltedStressPs = MeltedStressPs %>% select(Strain,GeneAdjust,sigLabel)

MeltedStressFC = MeltedStress %>% filter(variable=="log2FoldChange")
MeltedStressFC = MeltedStressFC %>% left_join(MeltedStressPs, by=c("Strain", "GeneAdjust"))

ggplot(MeltedStressFC, aes(x=Strain, y=GeneAdjust, fill=value)) + geom_tile() + scale_fill_gradient2(high="red", mid="white", low="blue", midpoint = 0) + geom_text(aes(label=sigLabel)) + theme_classic() + labs(fill="log2FC after H202") + ylab("Gene")


# Gene set enrichment analysis for GO terms
################################################

# Which gene ontologies are enriched in genes upregulated in DORN925 in response to stress? 
geneList
MakeGOMapForCP = genemap %>% select(GO_terms,AnnotID)
MakeGOMapForCP$GO_terms = sapply(MakeGOMapForCP$GO_terms, function(x) str_remove_all(x," "))
UncollapsedGOs = tidyr::separate_rows(MakeGOMapForCP,GO_terms, sep=";" )

UncollapsedGOs = data.frame(UncollapsedGOs)
UncollapsedGOs = UncollapsedGOs %>% filter(GO_terms!= "" & !is.na(GO_terms))
UncollapsedGOs


# DORN925 Stress to Non Stress 
enrichGO
go2ont(CPGo$GO)

termtogene = UncollapsedGOs

# Gene list for input into GSEA()/enricher()

MyGeneList = results_DORN925_StressDF$log2FoldChange
names(MyGeneList) <- results_DORN925_StressDF$AnnotID

MyGeneList1088 = results_DORN1088_StressDF$log2FoldChange
names(MyGeneList1088) = results_DORN1088_StressDF$AnnotID
MyGeneList1088 = rev(sort(MyGeneList1088))


MyGeneList = rev(sort(MyGeneList))

DORN925StressGSEA = GSEA(MyGeneList, TERM2GENE = termtogene)
DORN1088StressGSEA = GSEA(MyGeneList1088, TERM2GENE = termtogene)
DORN1088StressEnrichmentDown= enricher(DOWN_1088$AnnotID, TERM2GENE = termtogene)

enrichplot::ridgeplot(DORN1088StressGSEA) + ggtitle("GO Terms Up/Downregulated in DORN1088 in response to stress(GSEA)")
# UPREGULATED: GO:0006099 tricarboxylic acid cycle
# DOWNREGULATED:
# GO:0005886: plasma membrane
# GO:0043190 ATP-binding cassette (ABC) transporter complex
# GO:0005840 ribosome
# GO:0006189  'de novo' IMP biosynthetic process
# https://www.nature.com/articles/s41467-019-08724-x

enrichplot::ridgeplot(DORN925StressGSEA) + ggtitle("GO Terms Up/Downregulated in DORN925 in response to stress(GSEA)")

DORN925StressEnrichmentUp= enricher(UP_925$AnnotID, TERM2GENE = termtogene)
# "GO:0006099" "GO:0016829" "GO:0015074"
# tricarboxylic acid cycle
DORN925StressGSEA

DORN925StressEnrichmentDown= enricher(DOWN_925$AnnotID, TERM2GENE = termtogene)
View(DORN1088StressEnrichmentDown@result)

barplot(DORN925StressEnrichment)

enrichplot::heatplot(DORN925StressGSEA, foldChange=MyGeneList)

dotplot(DORN925StressGSEA, showCategory=30) + ggtitle("dotplot for GSEA")


# GO terms enriched in up-regulated genes in DORN925 in response to stress :
#############################################################################
(results_DORN925_StressDF %>% filter(grepl(GO_terms, pattern="GO:0006099")))$GeneAdjust
# Genes with the GO:0006099 term:
# "fumC", "Q614_SASC00275G0001", "SAOUHSC_01801", "odhA",
# "odhB", "SAOUHSC_01347", "sucD", "sucC", "SAOUHSC_01105", "mqo", "mqo" 


(results_DORN925_StressDF %>% filter(grepl(GO_terms, pattern="GO:0016829")))$GeneAdjust
# Genes in GO:0016829
#"SAOUHSC_01870", "hemB", "SAOUHSC_01587", "SAOUHSC_01307"(ltaE), "purE", "eno"
# "SAOUHSC_00705", "hisH", "hisF", "SAOUHSC_02685", "fba"


(results_DORN925_StressDF %>% filter(grepl(GO_terms, pattern="GO:0015074")))$GeneAdjust





DEstrain_control = results_genotypeDF %>% filter(padj <= .05)
DEstrain_stress = results_genotypeH202DF %>% filter(padj <= .05)


DEstrain_control$Direction = if_else(DEstrain_control$log2FoldChange<0, "Down", "Up")
DEstrain_control$DirectionAnnot = paste(DEstrain_control$AnnotID, DEstrain_control$Direction, sep="_")
DEstrain_controlTableDownTop20 = (DEstrain_control %>% filter(Direction=="Down") %>% arrange(padj) %>% select(padj, log2FoldChange, Gene))[1:20,]
DEstrain_controlTableUpTop20 = (DEstrain_control %>% filter(Direction=="Up") %>% arrange(padj) %>% select(padj, log2FoldChange, GeneAdjust))[1:20,]

write.csv(DEstrain_controlTableDownTop20, file='/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/DownRegulatedControlTop20.csv')
write.csv(DEstrain_controlTableUpTop20, file='/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/UpRegulatedControlTop20.csv')

DEstrain_stress$Direction = if_else(DEstrain_stress$log2FoldChange<0, "Down", "Up")
DEstrain_stress$DirectionAnnot = paste(DEstrain_stress$AnnotID, DEstrain_stress$Direction, sep="_")
DEstrain_stressDownTop20 = (DEstrain_stress %>% filter(Direction=="Down") %>% arrange(padj) %>% select(padj, log2FoldChange, GeneAdjust))[1:20,]
DEstrain_stressUpTop20 = (DEstrain_stress %>% filter(Direction=="Up") %>% arrange(padj) %>% select(padj, log2FoldChange, GeneAdjust))[1:20,]


(DEstrain_stress %>% filter(Direction=="Down") %>% arrange(padj) %>% select(padj, log2FoldChange, GeneAdjust))[1:20,]
write.csv(DEstrain_stressDownTop20, file='/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/DownRegulatedStressTop20.csv')
DEstrain_stressUpTop20
write.csv(DEstrain_stressUpTop20, file='/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/UpRegulatedStressTop20.csv')

SharedDE_strain_BetweenConditions = intersect(DEstrain_control$AnnotID, DEstrain_stress$AnnotID )
SharedDE_strain_BetweenConditionsDirection = intersect(DEstrain_control$DirectionAnnot, DEstrain_stress$DirectionAnnot )

justAnnot = sapply(SharedDE_strain_BetweenConditionsDirection, function(x) paste((str_split(x,pattern="_"))[[1]][1], (str_split(x,pattern="_"))[[1]][2], sep="_") )
setdiff(SharedDE_strain_BetweenConditions, justAnnot)
setdiff(  SharedDE_strain_BetweenConditions,SharedDE_strain_BetweenConditionsDirection)

DEstrain_stressSubsetShared = DEstrain_stress %>% filter(DirectionAnnot %in% SharedDE_strain_BetweenConditionsDirection)
DEstrain_stressSubsetSharedJustAnnot = DEstrain_stressSubsetShared %>% filter(AnnotID %in% SharedDE_strain_BetweenConditions)


JustStress = setdiff( DEstrain_stress$AnnotID, DEstrain_control$AnnotID)
Justcontrol =  setdiff( DEstrain_control$AnnotID,DEstrain_stress$AnnotID)

SharedDE_control = DEstrain_control %>% filter(DirectionAnnot %in% intersect(DEstrain_control$DirectionAnnot, DEstrain_stress$DirectionAnnot ))
SharedDE = SharedDE_control %>% select(padj, log2FoldChange, AnnotID) %>% left_join(DEstrain_control,by="AnnotID")

write.csv(SharedDE, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/DEGenes_by_Strain_BothConditions.csv")


DEstrain_controlUp = DEstrain_control %>% filter(Direction=="Up")
DEstrain_controlDown = DEstrain_control %>% filter(Direction=="Down")

DEstrain_stressUp =  DEstrain_stress %>% filter(Direction=="Up")
DEstrain_stressDown =  DEstrain_stress %>% filter(Direction=="Down")


VennDiagramlist_Up = list(Control=DEstrain_controlUp$GeneAdjust, H202=DEstrain_stressUp$GeneAdjust)
names(VennDiagramlist_Up) = c("Control Conditions", "H202 Stress")
ggvenn::ggvenn(VennDiagramlist_Up, fill_color = c("#F59999","#F59999"))

VennDiagramlist_Down = list(Control=DEstrain_controlDown$GeneAdjust, H202=DEstrain_stressDown$GeneAdjust)
names(VennDiagramlist_Down) = c("Control Conditions", "H202 Stress")
ggvenn::ggvenn(VennDiagramlist_Down, fill_color = c("#B1D4E0","#B1D4E0"))





#######################################
# (iii) EXAMINING DE IN JUST CRT OPERON 
#######################################
# In just control conditions
allannotated = results_genotypeDF %>% select(padj,pvalue, log2FoldChange)
allannotated$AnnotID = rownames(allannotated)
allannotated= allannotated %>% left_join(MapGenes,by="AnnotID")
View(allannotated)
allannotated %>% filter(grepl(Gene, pattern="crt"))


# What about when pooled across both conditions?  
DEDataSTRAIN =DESeqDataSetFromMatrix(countData=CountsMatrix, 
                                     colData=columndataDE,
                                     design= ~Strain)

results_genotypeCrossCondition <- results(myDESeqObjSTRAIN, contrast=c("Strain", "DORN925","DORN1088"))
results_genotypeCrossConditionDF = data.frame(results_genotypeCrossCondition) %>% select(padj,pvalue, log2FoldChange)
results_genotypeCrossConditionDF$AnnotID = rownames(results_genotypeCrossConditionDF)
results_genotypeCrossConditionDF= results_genotypeCrossConditionDF %>% left_join(MapGenes,by="AnnotID")
View(results_genotypeCrossConditionDF)
results_genotypeCrossConditionDF %>% filter(grepl(Gene, pattern="crt"))

# What about in just the H202 conditions?
results_genotypeh202 <- results(myDESeqObj, contrast=c("group", "treatment_DORN925","treatment_DORN1088"))
results_genotypeh202DF = data.frame(results_genotypeh202) %>% select(padj,pvalue, log2FoldChange)
results_genotypeh202DF$AnnotID = rownames(results_genotypeh202DF)
results_genotypeh202DF= results_genotypeh202DF %>% left_join(MapGenes,by="AnnotID")
results_genotypeh202DF %>% filter(grepl(Gene, pattern="crt"))






