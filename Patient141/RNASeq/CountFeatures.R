# Amy Campbell
# Output of the BT alignment 

library(dplyr)
library(Rsubread)
library(Rsamtools)
library(stringr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)

BiocManager::install("clusterProfiler")
samplemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/MappingRNASeq.txt", sep="\t", header=F)
colnames(samplemap) = c("Sample", "SampleTreatment", "Condition")
genemap=read.csv2("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/GeneNamesPreliminary.csv", sep='\t')

genemap=read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/UniProtInfoMapIDs.tsv",sep='\t')
# bamdirpath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/bams/"
# bamdirectory = list.files(bamdirpath)
# bamdirectory = sapply(bamdirectory, function(x) paste0(bamdirpath, x))
# bamfiles=BamFileList(bamdirectory)
# seqinfo(bamfiles[1])
gtffile <- file.path("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/DORN925_ref.gff")

MyCounts = featureCounts(files = bamdirectory,annot.ext =gtffile, isGTFAnnotationFile=T, GTF.featureType="CDS", GTF.attrType = "ID", isPairedEnd=T)
MyCounts
save(MyCounts, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RNACounts.rda")

