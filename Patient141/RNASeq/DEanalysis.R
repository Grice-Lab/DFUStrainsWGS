library(dplyr)
library(Rsubread)
library(Rsamtools)
library(stringr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)

samplemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/MappingRNASeq.txt", sep="\t", header=F)
colnames(samplemap) = c("Sample", "SampleTreatment", "Condition")
genemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/GeneNamesPreliminary.csv", sep='\t')
load("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNACounts.rda")



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

# group information includes info about both condition & strain
columndataDE$group = factor(paste(columndataDE$Condition, columndataDE$Strain, sep="_"))

DEData =DESeqDataSetFromMatrix(countData=CountsMatrix, 
                               colData=columndataDE,
                               design= ~group)
myDESeqObj = DESeq(DEData)

###########################################################
# (i) Comparing DORN925 to DORN1088 in control conditions 
###########################################################
results_genotype <- results(myDESeqObj, contrast=c("group", "ctrl_DORN925","ctrl_DORN1088"))

MapGenes = data.frame(AnnotID=row.names(results_genotype))
MapGenes = MapGenes %>% left_join(genemap, by="AnnotID") %>% unique()

MAplotStrain = ggmaplot(results_genotype,size=1, genenames = MapGenes$Gene)
MAplotStrain = MAplotStrain + ggtitle("Differentially expressed genes between \nDORN925:DORN1088 under control conditions") + theme(plot.title=element_text(face="bold", size=20, hjust=.5))
ggsave(MAplotStrain, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotControlStrains.png")

results_genotypeDF = data.frame(results_genotype) %>% select(padj,pvalue, log2FoldChange)
results_genotypeDF$AnnotID = row.names(results_genotypeDF)
results_genotypeDF = results_genotypeDF %>% left_join(MapGenes, by="AnnotID")
results_genotypeDF %>% filter(AnnotID=="cds-pgaptmp_002685")


###########################################################
# (ii) Comparing DORN925 to DORN1088 in H202 conditions 
###########################################################
results_genotypeH202 <- results(myDESeqObj, contrast=c("group", "treatment_DORN925","treatment_DORN1088"))
results_genotypeH202DF = data.frame(results_genotypeH202) %>% select(padj,pvalue, log2FoldChange)

MapGenesH202 = data.frame(AnnotID = row.names(results_genotypeH202DF))
MapGenesH202 = MapGenesH202 %>% left_join(genemap, by="AnnotID") %>% unique()

MAplotStrainH202 = ggmaplot(results_genotypeH202,size=1, genenames = MapGenesH202$Gene) +  ggtitle("Differentially expressed genes between \nDORN925:DORN1088 after H202 exposure") + theme(plot.title=element_text(face="bold", size=20, hjust=.5))
ggsave(MAplotStrainH202, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNASeqPlots/MAPlotH202Strains.png")

results_genotypeH202DF$AnnotID = row.names(results_genotypeH202DF)
results_genotypeH202DF = results_genotypeH202DF %>% left_join(MapGenesH202, by="AnnotID")
results_genotypeCrossConditionDF= results_genotypeCrossConditionDF %>% left_join(MapGenes,by="AnnotID")






#################################
# EXAMINING DE IN JUST CRT OPERON 
#################################
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

myDESeqObjSTRAIN = DESeq(DEDataSTRAIN)
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



# Multifactor
####################
