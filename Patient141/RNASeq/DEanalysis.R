
# Amy Campbell
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


samplemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/MappingRNASeq.txt", sep="\t", header=F)
colnames(samplemap) = c("Sample", "SampleTreatment", "Condition")
genemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/GeneNamesPreliminary.csv", sep=',')
load("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNACounts.rda")

genemap = genemap %>% mutate(GeneAdjust= if_else( (is.na(Gene) | Gene==""), RefSeqID, Gene))

# Manipulate sample ID mapping file and make plot of Sample vs. Total Reads assigned
#####################################################################################
CountsBySample=data.frame((colSums(MyCounts$counts)))
CountsBySample$Sample=sapply(row.names(CountsBySample), function(x) str_split(x, pattern="_")[[1]][1])
colnames(CountsBySample) = c("TotalAssignedReads", "Sample")

CountsBySample = samplemap %>% left_join(CountsBySample, by="Sample")
CountsBySample[6,4 ] <- 0
CountsBySample$Strain = sapply(CountsBySample$SampleTreatment, function(x) if_else(grepl(x, pattern="ARM72"),"DORN925(high)","DORN1088(low)") )

CountsBySample$ConditionStrain = paste(CountsBySample$Strain, CountsBySample$Condition, sep=" : ")
CountsBySample$ConditionStrain = factor(CountsBySample$ConditionStrain )
TotalCountsPlot = ggplot(CountsBySample, aes(x=Sample, y=TotalAssignedReads, fill=ConditionStrain)) + geom_bar(stat="identity", color="black")
TotalCountsPlot$data$Sample = factor(TotalCountsPlot$data$Sample , levels=(TotalCountsPlot$data %>% arrange(ConditionStrain))$Sample)

TotalCountsPlot$data$Sample = factor(TotalCountsPlot$data$Sample , levels=(TotalCountsPlot$data %>% arrange(Strain, Condition))$Sample)
max(TotalCountsPlot$data$TotalAssignedReads)
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

columndataDE$group = factor(paste(columndataDE$Condition, columndataDE$Strain, sep="_"))

DEData =DESeqDataSetFromMatrix(countData=CountsMatrix, 
                               colData=columndataDE,
                               design= ~group)
myDESeqObj = DESeq(DEData)

# Sanity check:sak should be 0 counts in DORN1088 samples
##########################################################
CountsMatrix[c("cds-pgaptmp_000025"),]


##########################################
# How does DORN925 respond to H202 stress?
###########################################
results_DORN925_Stress <- results(myDESeqObj, contrast=c("group", "treatment_DORN925","ctrl_DORN925"))
results_DORN925_StressDF = (data.frame(results_DORN925_Stress))
results_DORN925_StressDF$AnnotID = row.names(results_DORN925_StressDF)
results_DORN925_StressDF = results_DORN925_StressDF %>% left_join(genemap, by="AnnotID")

MapGenes925 = data.frame(AnnotID = row.names(results_DORN925_StressDF))
MapGenes925 = MapGenes925 %>% left_join(genemap, by="AnnotID") %>% unique()

MAplot925 = ggmaplot(results_DORN925_StressDF,size=1, genenames = MapGenes925$GeneAdjust) + ggtitle("Differentially expressed genes in DORN925\nin Response to H202 Stress")
MAplot925= MAplot925 + theme(plot.title=element_text(face="bold", size=20, hjust=.5))

ggsave(MAplot925, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotDORN925Stress.png", height=)


crtGenes = c("cds-pgaptmp_002177", "cds-pgaptmp_002178", "cds-pgaptmp_002179", "cds-pgaptmp_002180", "cds-pgaptmp_002181")

CountsMatrix[crtGenes, ]


#########################################################
# DE in response to stress in DORN1088 (lowxanthin)
##########################################################

results_DORN1088_Stress <- results(myDESeqObj, contrast=c("group", "treatment_DORN1088","ctrl_DORN1088"))
results_DORN1088_StressDF = (data.frame(results_DORN1088_Stress))

MapGenes1088 = data.frame(AnnotID = row.names(results_DORN1088_StressDF))
MapGenes1088 = MapGenes1088 %>% left_join(genemap, by="AnnotID") %>% unique()

MAplot1088 = ggmaplot(results_DORN1088_StressDF,size=1, genenames = MapGenes1088$GeneAdjust) + ggtitle("Differentially expressed genes in DORN1088\nin Response to H202 Stress")

MAplot1088= MAplot1088 + theme(plot.title=element_text(face="bold", size=20, hjust=.5))

ggsave(MAplot1088, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotDORN1088Stress.png", height=)


# Volcano plot for 925
#########################
# Make log2fold change cutoff .6, representing 1.5X-fold change 
results_DORN925_StressDF$DE = if_else( abs(results_DORN925_StressDF$log2FoldChange) < .6  | results_DORN925_StressDF$padj >= .05 | is.na(results_DORN925_StressDF$padj), FALSE, TRUE)

results_DORN925_StressDF = results_DORN925_StressDF %>% mutate(DirectionDE = case_when( log2FoldChange > 0 & DE==TRUE ~ "Up",
                                                                                        log2FoldChange <0 & DE== TRUE ~ "Down", 
                                                                                        TRUE ~ "None"))
# only label those with 4X and higher (log2Foldchange>=2)
results_DORN925_StressDF = results_DORN925_StressDF %>% mutate(Label = if_else( ((DirectionDE == "Up" | DirectionDE == "Down" ) & abs(log2FoldChange) >2 ), GeneAdjust,""))
results_DORN925_StressDF$Label[results_DORN925_StressDF$Label==""] <- NA

VolcPlot925 = ggplot(data=results_DORN925_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=Label)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel() + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN925\nin Response to H202 Stress") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5))


# Volcano plot for 1088
#########################
results_DORN1088_StressDF$AnnotID = row.names(results_DORN1088_StressDF)
results_DORN1088_StressDF = results_DORN1088_StressDF %>% left_join(MapGenes1088, by="AnnotID")

# Make log2fold change cutoff .6, representing 1.5X-fold change 
results_DORN1088_StressDF$DE = if_else( abs(results_DORN1088_StressDF$log2FoldChange) < .6  | results_DORN1088_StressDF$padj >= .05 | is.na(results_DORN1088_StressDF$padj), FALSE, TRUE)

results_DORN1088_StressDF = results_DORN1088_StressDF %>% mutate(DirectionDE = case_when( log2FoldChange > 0 & DE==TRUE ~ "Up",
                                                                                          log2FoldChange <0 & DE== TRUE ~ "Down", 
                                                                                          TRUE ~ "None"))
# only label those with 4X and higher (log2Foldchange>=2)
results_DORN1088_StressDF = results_DORN1088_StressDF %>% mutate(Label = if_else( ((DirectionDE == "Up" | DirectionDE == "Down" ) & abs(log2FoldChange) >2 ), GeneAdjust, ""))
results_DORN1088_StressDF$Label[results_DORN1088_StressDF$Label==""] <- NA


VolcPlot1088 = ggplot(data=results_DORN1088_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=Label)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel() + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN1088\nin Response to H202 Stress") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5))





gridarrange= gridExtra::grid.arrange(VolcPlot1088, VolcPlot925, ncol=2)
ggsave(gridarrange, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/SideBySide.pdf", width=20, height=10)


UP_1088 = results_DORN1088_StressDF %>% filter(DirectionDE=="Up")
DOWN_1088 = results_DORN1088_StressDF %>% filter(DirectionDE=="Down")

UP_925 = results_DORN925_StressDF %>% filter(DirectionDE=="Up")
DOWN_925 = results_DORN925_StressDF %>% filter(DirectionDE=="Down")

UniqueUp925 = results_DORN925_StressDF %>% filter(AnnotID %in% setdiff(UP_925$AnnotID, UP_1088$AnnotID))
UniqueDown925 = results_DORN925_StressDF %>% filter( (AnnotID %in% setdiff(DOWN_925$AnnotID, DOWN_1088$AnnotID)) )


UniqueUp1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% setdiff(UP_1088$AnnotID,UP_925$AnnotID) )
UniqueDown1088 = results_DORN1088_StressDF %>% filter(AnnotID %in% setdiff(DOWN_1088$AnnotID,DOWN_925$AnnotID) )



results_DORN925_StressDF = results_DORN925_StressDF %>% mutate(LabelUnique = if_else(DE==TRUE & (AnnotID %in%   append(UniqueUp925$AnnotID, UniqueDown925$AnnotID)),GeneAdjust, ""))
results_DORN925_StressDF$LabelUnique[results_DORN925_StressDF$LabelUnique ==""] <- NA
VolcPlot925UniqueGr2 = ggplot(data=results_DORN925_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=LabelUnique)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel(color="black") + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN925\nin Response to H202 Stress (Unique)") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5))



results_DORN1088_StressDF = results_DORN1088_StressDF %>% mutate(LabelUnique = if_else(DE==TRUE & (AnnotID %in%   append(UniqueUp1088$AnnotID, UniqueDown1088$AnnotID)),GeneAdjust, ""))
results_DORN1088_StressDF$LabelUnique[results_DORN1088_StressDF$LabelUnique ==""] <- NA
VolcPlot1088UniqueGr2 = ggplot(data=results_DORN1088_StressDF, aes(x=log2FoldChange, y=-log10(pvalue), col=DirectionDE, label=LabelUnique)) + 
  geom_point() + 
  theme_minimal() +geom_text_repel(color="black") + scale_color_manual(values=c("royalblue4", "black", "firebrick3")) +  geom_vline(xintercept=c(-0.6, 0.6), col="grey52",linetype="dashed") +geom_vline(xintercept=c(-2, 2), col="grey52",linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey52", linetype="dashed") + geom_hline(yintercept=-log10(0.01), col="grey52", linetype="dashed")  + ggtitle("Differentially expressed genes in DORN1088\nin Response to H202 Stress(Unique)") + theme_minimal() + theme(plot.title=element_text(face="bold", size=20, hjust=.5))


arranged = gridExtra::grid.arrange(VolcPlot1088UniqueGr2,VolcPlot925UniqueGr2,ncol=2)
ggsave(arranged, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/SideBySideUnique.pdf", width=20, height=10)
###########################################################
# (i) Comparing DORN925 to DORN1088 in control conditions 
###########################################################
# "group" is the combination of strain and condition 


results_genotype <- results(myDESeqObj, contrast=c("group", "ctrl_DORN925","ctrl_DORN1088"))

MapGenes = data.frame(AnnotID=row.names(results_genotype))
MapGenes = MapGenes %>% left_join(genemap, by="AnnotID") %>% unique()

MAplotStrain = ggmaplot(results_genotype,size=1, genenames = MapGenes$GeneAdjust)
MAplotStrain = MAplotStrain + ggtitle("Differentially expressed genes between \nDORN925:DORN1088 under control conditions") + theme(plot.title=element_text(face="bold", size=20, hjust=.5))
ggsave(MAplotStrain, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotControlStrains.png")

results_genotypeDF = data.frame(results_genotype) %>% select(padj,pvalue, log2FoldChange)
results_genotypeDF$AnnotID = row.names(results_genotypeDF)
results_genotypeDF = results_genotypeDF %>% left_join(MapGenes, by="AnnotID")


###########################################################
# (ii) Comparing DORN925 to DORN1088 in H202 conditions 
###########################################################
results_genotypeH202 <- results(myDESeqObj, contrast=c("group", "treatment_DORN925","treatment_DORN1088"))
results_genotypeH202DF = data.frame(results_genotypeH202) %>% select(padj,pvalue, log2FoldChange)

MapGenesH202 = data.frame(AnnotID = row.names(results_genotypeH202DF))
MapGenesH202 = MapGenesH202 %>% left_join(genemap, by="AnnotID") %>% unique()

MAplotStrainH202 = ggmaplot(results_genotypeH202,size=1, genenames = MapGenesH202$GeneAdjust) +
  ggtitle("Differentially expressed genes between \nDORN925:DORN1088 after H202 exposure") + theme(plot.title=element_text(face="bold", size=20, hjust=.5))
ggsave(MAplotStrainH202, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotH202Strains.png")

results_genotypeH202DF$AnnotID = row.names(results_genotypeH202DF)
results_genotypeH202DF = results_genotypeH202DF %>% left_join(MapGenesH202, by="AnnotID")




# Which genes are DE between low-high only in control? Only under stress? under both conditions?

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


# c("#B1D4E0", "#79A9F5")
s
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






