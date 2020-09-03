# Amy Campbell
# November 2019
# Shotgun Dirichlet Analysis
# Reproducing Loesche's analysis from 
# "Temporal Stability in Chronic Wound Microbiota Is Associated With Poor Healing"
# This time in the 46 Subjects for whom we have metagenomic shotgun data 
# Run originally from Club_Grice/ working directory 


# Experiment: 
##############
# "Are there distinct community types within the shotgun data for the DFU cohort, and are they the same as in the 16S?"

# Required software
###################
#BiocManager::install("DirichletMultinomial")
library(dplyr)
library(lattice)
library(tidyr)
library(ggplot2)
library(DirichletMultinomial)

# Color palettes
#################
library(awtools)
library(ggthemes)
library(stringi)

# Input Files
#############
AbundancePath = "scripts/acampbe/DFU/Abundance_FivePercentCutoff.csv"
Abundance = read.csv(AbundancePath, row.names=1)
ReadsPerSample = "scripts/acampbe/DFU/scripts/shotgun_analysis_scripts/ReadCountsTimepoints.csv"
BacterialPctSample = "scripts/acampbe/DFU/scripts/shotgun_analysis_scripts/BacteriaPct.csv"
LoescheResults = "scripts/loesche/dfu_uclust/dmn/dmn_results.csv"
LoescheMapping="scripts/loesche/dfu_uclust/metadata/v1v3_sample_key.csv"

# I'm going to read this in from the file output by TimeSeriesSummary.R


# These are data frames which say how many reads were in each fastq before input into MetaPhlan(e.g. after QC), 
# As well as the % of reads metaphlan predicted to be bacterial in that sample
Readcts = read.csv(ReadsPerSample, header=0)
Bacteria = read.csv(BacterialPctSample, header=0)
colnames(Bacteria) = c("SampleName", "PctBact")
colnames(Readcts) = c("SampleName", "Reads")

Readcts$SampleName <- sapply(Readcts$SampleName, function(x) stringi::stri_replace_last(str=x, replacement = '.', fixed="_"))
Bacteria$SampleName <- sapply(Bacteria$SampleName, function(x) stringi::stri_replace_last(str=x, replacement = '.', fixed="_"))

####
  
Abundance_Dirichlet = Abundance %>% select(-SubjectID, -Timepoint, -healed_by_12)


# Loesche previously filtered to things present in at least 10% of samples
# Though i'm always unsure whether to include things present at < .1% as 'present', he did so I guess I will. 
# Especially given that he was filtering on count data and I am now filtering on 'estimated count data'

Abundance_Dirichlet$SampleName = rownames(Abundance_Dirichlet)
Abundance_Dirichlet = Abundance_Dirichlet %>% left_join(Readcts, by="SampleName")
Abundance_Dirichlet =  Abundance_Dirichlet %>% left_join(Bacteria, by="SampleName")
Abundance_Dirichlet$PctBact = .01*Abundance_Dirichlet$PctBact
Abundance_Dirichlet <-(Abundance_Dirichlet[,(colSums(Abundance_Dirichlet >= 0) > 19)])
rownames(Abundance_Dirichlet) = Abundance_Dirichlet$SampleName
######

ReadNo = Abundance_Dirichlet$Reads*Abundance_Dirichlet$PctBact
Abundance_Dirichlet$PctBact = NULL
Abundance_Dirichlet$Reads = NULL
Abundance_Dirichlet$SampleName = NULL
Abundance_Dirichlet_Counts = apply(Abundance_Dirichlet, 2, function(x) round(ReadNo*x))

Abundance_Dirichlet
colSums(Abundance_Dirichlet_Counts)

# this leaves 49 taxa
dim((Abundance_Dirichlet_Counts))

counts = data.matrix((Abundance_Dirichlet_Counts))
row.names(counts) 
densityplot(log10(colSums(counts)), xlim=range(log10(colSums(counts))),xlab="Taxon representation (log 10 count)")
densityplot((colSums(counts)), xlim=range((colSums(counts))),xlab="Taxon representation (raw count)")
Abundance_Dirichlet_Counts_Save = Abundance_Dirichlet_Counts
# Fitting up to 8 components because I was burned whilst testing a wrong version of this model and was fooled by a local minimum 
FitDirichlet  = mclapply(1:8, dmn, count=counts, verbose=TRUE, seed=19143)
lplc <- sapply(FitDirichlet, laplace)

plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Model Fit', main="Dirichlet Components by Laplace Model Fit (Read Count Estimates)") 
ggplot2::qplot(x=c(1:8), y = lplc) + geom_line() + xlab("# Dirichlet Components Fit") + ylab("Model Fit(Laplac Approximation)") + ggtitle("Dirichlet Components by Laplace Model Fit (Read Count Estimates)")
best <- FitDirichlet[[which.min(lplc)]]
# K = 4 
df = data.frame(best@group)
colnames(df) = c("Group1", "Group2", "Group3", "Group4")

groups = data.frame(best@group)
colnames(groups)= c("Group1", "Group2", "Group3", "Group4")
groups$cluster <- sapply(1:nrow(groups), FUN=function(x) which.max(groups[x,]))

Abundance_Dirichlet_Counts <- data.frame(Abundance_Dirichlet_Counts)
groups$SampleID = rownames(groups)
Abundance_Dirichlet_Counts$SampleID = rownames(Abundance_Dirichlet_Counts)
Abundance_Dirichlet_Counts = Abundance_Dirichlet_Counts %>% left_join(groups[c("cluster", "SampleID")], by="SampleID")


Abundance_Dirichlet_Proportions = Abundance_Dirichlet_Counts
abundance_prop_intrmediate = Abundance_Dirichlet_Proportions[,1:(ncol(Abundance_Dirichlet_Proportions)-2)] 

Abundance_Dirichlet_Proportions[,1:(ncol(Abundance_Dirichlet_Proportions)-2)] =  t(apply(abundance_prop_intrmediate, 1, function(x) x/(sum(x))))


# Put healing info, SubjectID, Timepoint back into the dataframe
Abundance$SampleID = rownames(Abundance)

Abundance_Dirichlet_Proportions = Abundance_Dirichlet_Proportions %>% full_join(Abundance[c("SampleID","SubjectID","Timepoint", "healed_by_12")], by="SampleID")
means = Abundance_Dirichlet_Proportions %>% group_by(cluster) %>% select(-SampleID,-SubjectID,-Timepoint,  -healed_by_12 ) %>% summarize_all(mean)

means= data.frame(means)
clustersave = means$cluster
means$cluster = NULL
colnames(means)[apply(means[, 1:ncol(means)],1,which.max)]

Abundance_Dirichlet_Proportions
View(Abundance_Dirichlet_Proportions)
Save_Factors  = Abundance_Dirichlet_Proportions[c("SampleID", "cluster", "SubjectID", "Timepoint", "healed_by_12")]
Abundance_Dirichlet_Proportions = Abundance_Dirichlet_Proportions %>% select(-SampleID, -cluster, -SubjectID, -Timepoint, -healed_by_12)
#Abundance_Dirichlet_Proportions <- (Abundance_Dirichlet_Proportions[,(colSums(Abundance_Dirichlet_Proportions >= .05) > 10)])

Abundance_Dirichlet_Proportions = cbind(Save_Factors,Abundance_Dirichlet_Proportions )
abundance_melt = reshape2::melt(Abundance_Dirichlet_Proportions, id=c("SampleID","cluster", "SubjectID","Timepoint", "healed_by_12"))

subset  = abundance_melt %>% group_by(variable) %>%summarize(mean1=mean(value)) %>% arrange(-mean1)
subset_vars = subset[1:16,]$variable

top6 = abundance_melt %>% group_by(variable, cluster) %>% summarize(val =mean(value)) %>% group_by(cluster)%>% top_n(., 6, val)

unique(top6$variable)
abundance_melt = abundance_melt %>% filter(variable %in% top6$variable)
abundance_melt = abundance_melt %>% arrange(value)
#%>%  arrange(value) %>% top_n(.,10, value) %>% arrange((variable))
modified_simpsons = c("#FED439FF","#709AE1FF","#197EC0FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF","#9A009A","#71D0F5FF","#370335FF","#46732EFF", "#075149FF","#C80813FF","#91331FFF", "#1A9993FF", "#FD8CC1FF")

coveragestats = read.csv("/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/scripts/shotgun_analysis_scripts/CoverageStats.tsv", sep='\t')
coveragestats$PctCovered = coveragestats$PctCovered / 2782740

View(coveragestats)
healeddf %>% healedD

qplot(coveragestats$PctCovered, coveragestats$TotalDepth)
#ggplot() + geom_bar(data=abundance_melt, aes(x=cluster, y=value, fill=variable), stat="identity") + scale_fill_manual(values=modified_simpsons) + xlab("Cluster Identity") +
#  ylab("Sum Relative Abundance") +
#  ggtitle("Distribution of Top 6 Abundant Species by Dirichlet Component") + 
#  theme(legend.title = element_text(color = "black", size = 16), legend.text =element_text(size=12), title=element_text(size=20, face='bold'), axis.text.x = element_text(size=14, face='bold'), axis.text.y = element_text(size=12) ) + labs(fill="Species")

top6 = top6 %>% arrange(val)
ggplot() + geom_bar(data=top6, aes(x=cluster, y=val, fill=variable), stat="identity", position="stack") + scale_fill_manual(values=modified_simpsons) + xlab("Cluster Identity") +
  ylab("Average Relative Abundance") +
  ggtitle("Distribution of Top 6 Abundant Species by Dirichlet Component") + 
  theme(legend.title = element_text(color = "black", size = 16), legend.text =element_text(size=12), title=element_text(size=20, face='bold'), axis.text.x = element_text(size=14, face='bold'), axis.text.y = element_text(size=12) ) + labs(fill="Species")



abundance_melt$value
scale_fill_manual()
abundance_melt_updated = abundance_melt %>% group_by(variable, cluster) %>% summarize(val=mean(value))
# Simspons palette default(ggsci)
ggplot() + geom_bar(data=abundance_melt_updated, aes(x=cluster, y=val, fill=variable), stat="identity", position="stack") + scale_fill_manual(values=modified_simpsons) + ylab("Mean Relative Abundance") + xlab("Cluster Identity") + ggtitle("Most abundant species in dirichlet components from DFU shotgun data")


plot(Abundance_Dirichlet_Proportions$cluster, Abundance_Dirichlet_Proportions$healed_by_12)

Loesche = read.csv(LoescheResults)
Loeschemap = data.frame(read.csv(LoescheMapping))
View(Loeschemap)

Loeschemap$SubjectID <- sapply(Loeschemap$SubjectID, function(x) stringr::str_replace(x, "-", "."))
Loeschemap$SubjectID <- sapply(Loeschemap$SubjectID, function(x) paste0("filtered_sorted_", x))

Loeschemap$SampleID = as.factor(Loeschemap$SampleID)
Loesche$SampleID = as.factor(Loesche$SampleID)
Loeschemap = Loeschemap %>% full_join(Loesche, by="SampleID")

Loeschemap = Loeschemap[c("SubjectID", "dmn")]
colnames(Loeschemap) = c("SampleID", "Old_Cluster")
New_Proportions = Abundance_Dirichlet_Proportions %>% left_join(Loeschemap, by="SampleID")

View(New_Proportions)
New_Proportions = New_Proportions %>% filter(!is.na(Old_Cluster))
plot(New_Proportions$cluster, New_Proportions$Old_Cluster)
View(New_Proportions[c("SampleID","cluster", "Old_Cluster")])
ggplot() + geom_bar(data=New_Proportions, aes(x=Old_Cluster, fill=cl), stat="prop", group="Old_Cluster") + scale_fill_manual(values=modified_simpsons)
ggplot(New_Proportions, aes(x=Old_Cluster, group=as.factor(cluster))) + geom_bar(aes(fill=as.factor(cluster)), position="fill") +
  scale_fill_manual(values=rev(awtools::spalette)) + ggtitle("Overlap Between 16S DMM clusters and Shotgun DMM clusters in DFU") + xlab("16S-based component") + ylab("Proportion") + labs(fill="Shotgun-based component")

ggplot(Abundance_Dirichlet_Proportions, aes(x=as.factor(SubjectID), group=as.factor(cluster))) + geom_bar(aes(fill=as.factor(cluster)), position="fill") + scale_fill_manual(values=rev(modified_simpsons)) +facet_grid(.~healed_by_12)+ ggtitle("Proportions of 16S-based Dirichlet Assignments in Shotgun-based ")
ordered = (Abundance_Dirichlet_Proportions %>% group_by(SubjectID) %>% arrange(healed_by_12))$SampleID

timeplot = ggplot2::ggplot(Abundance_Dirichlet_Proportions, aes(y=factor(SubjectID), x=as.factor(Timepoint), fill = factor(cluster))) + geom_tile(color="white") + scale_fill_manual(values =rev(awtools::spalette))+ggtitle("Dirichlet community type by subject over time") + labs(fill="Cluster") + ylab("Subject ID") + theme_clean() + theme(plot.title=element_text(size=20, face='bold'), axis.title.x = element_text(face='bold', size=15), axis.title.y = element_text(face='bold', size=15), axis.text.x = element_text(size=14),axis.text.y = element_text(size=14), legend.position="bottom", legend.title=element_text(size=15),legend.text=element_text(size=14)) + xlab("Timepoint (separated by 2 weeks)")

timeplot = ggplot2::ggplot(Abundance_Dirichlet_Proportions, aes(y=factor(SubjectID), x=as.factor(Timepoint), fill = factor(cluster))) + geom_tile(color="white") + scale_fill_manual(values =rev(awtools::spalette))+ggtitle("Dirichlet community type by subject over time") + labs(fill="Cluster") + ylab("Subject ID") + theme_clean() + theme(plot.title=element_text(size=19, face='bold'), axis.title.x = element_text(face='bold', size=15), axis.title.y = element_text(face='bold', size=15), axis.text.x = element_text(size=14),axis.text.y = element_text(size=14), legend.position="bottom", legend.title=element_text(size=15),legend.text=element_text(size=14)) + xlab("Timepoint (separated by 2 weeks)")

timeplot2 = ggplot2::ggplot(New_Proportions, aes(y=factor(SubjectID), x=as.factor(Timepoint), fill = factor(cluster))) + geom_tile(color="white") + scale_fill_manual(values =rev(awtools::spalette))+ labs(fill="Cluster") + ylab("Subject ID") + theme_clean() + theme(plot.title=element_text(size=20, face='bold'), axis.title.x = element_text(face='bold', size=15), axis.title.y = element_text(face='bold', size=15), axis.text.x = element_text(size=14),axis.text.y = element_text(size=14), legend.position="bottom", legend.title=element_text(size=15),legend.text=element_text(size=14)) + ggtitle("Shotgun-Based")+ xlab("Timepoint (separated by 2 weeks)")
timeplot2_old = ggplot2::ggplot(New_Proportions, aes(y=factor(SubjectID), x=as.factor(Timepoint), fill = factor(Old_Cluster))) + geom_tile(color="white") + scale_fill_manual(values = (awtools::mpalette))+ labs(fill="Cluster") + ylab("Subject ID") + theme_clean() + theme(plot.title=element_text(size=20, face='bold'), axis.title.x = element_text(face='bold', size=15), axis.title.y = element_text(face='bold', size=15), axis.text.x = element_text(size=14),axis.text.y = element_text(size=14), legend.position="bottom", legend.title=element_text(size=15),legend.text=element_text(size=14))+ggtitle("16S-Based") + xlab("Timepoint (separated by 2 weeks)")
gridExtra::grid.arrange(timeplot2, timeplot2_old, ncol=2, top="Dirichlet components over time 16S vs. shotgun-based DFU metagenomes")

healed = Abundance_Dirichlet_Proportions %>% filter(healed_by_12==0) %>% arrange(SubjectID,Timepoint)
unhealed =Abundance_Dirichlet_Proportions %>% filter(healed_by_12==1) %>% arrange(SubjectID, Timepoint)
healed$SampleID
order = c(unique(healed$SubjectID), unique(unhealed$SubjectID))
timeplot
timeplot$data$SubjectID = factor(timeplot$data$SubjectID, levels=order)

ggsave("TestingDirichlet.png" , plot = timeplot, dpi=500, height=12, width=7)


for_diversity = Abundance_Dirichlet_Proportions %>% select(-SampleID, -SubjectID, -Timepoint, -healed_by_12)

for_diversity %>% group_by(cluster) %>% summarize_each( vegan::diversity(., index="shannon"))
Abundance_Dirichlet_Proportions$alpha =  vegan::diversity(for_diversity, index="shannon")

ggplot(Abundance_Dirichlet_Proportions, aes(x=as.factor(cluster), y=alpha, fill=as.factor(cluster))) + geom_boxplot() + 
  scale_fill_manual(values =rev(awtools::spalette)) + xlab("Dirichlet Cluster") + ylab("Shannon Diversity") + ggtitle("Shannon Diversity in Dirichlet Community Members") + 
  theme(axis.title.x=element_text(size=16, face='bold'), axis.title.y=element_text(size=16, face='bold'),axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),legend.title=element_text(size=15),legend.text=element_text(size=14), plot.title=element_text(size=20, face='bold'), legend.position="none") 
attach(Abundance_Dirichlet_Proportions)
pairwise.t.test(alpha, cluster, p.adj="bonf")
detach(Abundance_Dirichlet_Proportions)
