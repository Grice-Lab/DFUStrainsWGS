# Amy Campbell
# 05/2020

# Subsetting the DOERN S. aureus isolates to 'likely sisters' 
# for the purposes of doing lower-throughput phenotypic analyses, planning
# long read sequencing
# Basically, I'll need to come up with 'clusters' of likely sister isolates. 
# From there, we'll need to decide which one to use as the 'representative' one (perhaps based on length of total cleaned genome? coverage of it? )

# Read in pair-wise SNP distance between core genome alignments (alignments generated in MAFFT via roary, 
# SNP distances calculated using snp-dists from the makers of Prokka) 


#library("viridis")
library("ggplot2")
library("reshape2")
library("dplyr")

# Read in TSV of SNP distances 
SNPdist = read.csv("SNPDistances_DORNs.tsv", sep="\t", header=1, row.names=1)

#Mapping file for the DORN isolates to visit and patient ID
StaphIsolateDORNs=read.csv("DFU_Staph_aureus_isolates.csv")

rownames(SNPdist) = colnames(SNPdist)

#StaphIsolateDORNs$visit_1 = StaphIsolateDORNs$visit 
StaphIsolateDORNs$subject_timepoint = paste(StaphIsolateDORNs$patient_id, StaphIsolateDORNs$visit,sep="_")

StaphIsolateDORNs$Genome2=paste0("DORN", StaphIsolateDORNs$Doern.lab.bank.)
StaphIsolateDORNs$Genome1=paste0("DORN", StaphIsolateDORNs$Doern.lab.bank.)

# Make one column to add the 'genome 1' associated visit
# Make another to add the 'genome 2' associated visit 
for_merge = StaphIsolateDORNs[c("Genome1","patient_id", "visit", "subject_timepoint")]
for_merge2 = StaphIsolateDORNs[c("Genome2","patient_id", "visit")]
# make sure these new columns are clearly tied to the Genome 2 
colnames(for_merge2) = c("Genome2", "patient_id_2", "visit2")

SNPdist$firstgenome=row.names(SNPdist)
Melted = reshape2::melt(SNPdist, id.vars=c("firstgenome"))
colnames(Melted) = c("Genome1", "Genome2", "Distance")


Melted = Melted %>% left_join(for_merge, by="Genome1")
Melted = Melted %>% left_join(for_merge2, by="Genome2")
Melted = Melted %>% mutate(Almost_Identical = (Distance < 118) )
Melted = Melted %>% mutate(Identical = (Distance==0) )
table(Melted$patient_id)

# Assign Timepoints to Viridis "Magma" colors so we can code the 'ticks' colors by them (so we only have to explicitly label subject on the axes)
# I never actually use this so may get rid of it, but in case I do for some reason decide to do it again
FourteenTimepoints=viridis::magma(14)
visit_color= Melted$visit
for(colornum in 0:13){
  visit_color[which(visit_color==colornum)] <- FourteenTimepoints[(colornum+1)]
}

# Labeling 
Melted$visit_color = visit_color


# Plot 1  : "Almost identical" core genomes (>99.99% SNP identity)
#####################################################################
p1 <- ggplot(data = Melted, aes(x=Genome2, y=Genome1, fill=Almost_Identical)) + 
  geom_tile(color = "Black") + labs(fill=">99.99% single \ncore nucleotide identity") + 
  scale_fill_manual(values=c("#0072B2", "#F0D642")) +
  ggtitle("Binary Matrix Showing Isolates with <.01% SNP Distance Between Core Genomes") + 
  theme( plot.title=element_text(hjust=.5), axis.text.x=element_text(size=7, angle=90))

# Labels and reordering of plot p1
##################################
# Order rows of a DF representing the data in plot p1 by subject_timepoint
plotDF = data.frame(p1$data)
LevelsGenome1 = unique((plotDF %>% arrange(subject_timepoint))$Genome1)

# Hacky way to get ordering of subject timepoint associated with this ordering 
# After ordering the Genome1 names by their associated subject_timepoint order, 
# Build a list of labels for the axes corresponding to the correct subject_timepoint values
# This has to be done somewhat manually unfortunately, but that's ok.. 

# Make the DF
subjectTimepointDF = data.frame(LevelsGenome1)
colnames(subjectTimepointDF) = c("Genome1")
MappingSubjectTimepoint_order = subjectTimepointDF %>% left_join(StaphIsolateDORNs[c("Genome1", "subject_timepoint", "patient_id")])
Subj_tp_order = MappingSubjectTimepoint_order$subject_timepoint
subjORder=MappingSubjectTimepoint_order$patient_id

Greys=c("#D0D0D0", "#707070")
Grey_Assignments = rep(Greys, ceiling(length(unique(MappingSubjectTimepoint_order$patient_id))/2 ))
Subjects_Unique = unique(MappingSubjectTimepoint_order$patient_id)
SubjectColor=MappingSubjectTimepoint_order$patient_id

for(i in 1:length(Subjects_Unique)){
  SubjectColor[which(SubjectColor==Subjects_Unique[i])] <- Grey_Assignments[i]
}
p1$data$Genome2 = factor(p1$data$Genome2, levels=LevelsGenome1)
p1$data$Genome1 = factor(p1$data$Genome1, levels=LevelsGenome1) # 

# Save plot p1 as a file
########################
postscript(file="SimilarIdentityMatrix.eps", width=5000, height=5000)
p1 + geom_tile(color = "black", size=.001) +
  theme(axis.ticks.y = element_line(size=1.1, color=SubjectColor),
        axis.ticks.x = element_line(size=1.1, color=SubjectColor),
        axis.ticks.length.x = unit(1, "cm"),
        axis.ticks.length.y = unit(1, "cm"),
        plot.title=element_text(hjust=.5),
        axis.text.x=element_text(size=3, angle=90, hjust=1),
        axis.text.y=element_text(size=3, vjust=.5)) +
  scale_x_discrete(labels=subjORder) + scale_y_discrete(labels=subjORder) + coord_fixed(ratio=1)

dev.off()  

# Make plot p2 using the reordered levels from p1 but instead of filling by 'AlmostIdentical' factor (<118 SNPs diff),
# Fill by exact identical (0 SNP distance) similarity 
#################################################################################################################
p2 <- ggplot(p1$data, aes(x=Genome1, y=Genome2, fill=Identical)) + geom_tile(color = "black", size=.001) +
  theme(axis.ticks.y = element_line(size=1.1, color=SubjectColor),
        axis.ticks.x = element_line(size=1.1, color=SubjectColor),
        axis.ticks.length.x = unit(1, "cm"),
        axis.ticks.length.y = unit(1, "cm"),
        plot.title=element_text(hjust=.5),
        axis.text.x=element_text(size=3, angle=90, hjust=1),
        axis.text.y=element_text(size=3, vjust=.5)) +
  scale_fill_manual(values=c("#0072B2", "#F0D642")) +
  labs(fill=">Identical single \ncore nucleotide identity") + 
  scale_x_discrete(labels=subjORder) + scale_y_discrete(labels=subjORder) + coord_fixed(ratio=1) + 
  ggtitle("Binary Matrix Showing Isolates with 0 SNP Distance Between Core Genomes")

postscript(file="IdentityMatrixIdentical.eps", width=5000, height=5000)
p2
dev.off()

# Take away the rows that represent one isolate's match to itself 
Melted = Melted %>% mutate(LessThan5_Snps = (Distance < 5))
Melted_DifferentGenomePairs = Melted %>% filter(!(Genome1==Genome2))

# Because we don't want to, for the purposes of Oxford long read sequencing, 
# consider any isolates from different subjects 
# to be 'redundant,' we'll now only cluster within subjects 
Melted_DifferentGenomePairs_SameSubj = Melted_DifferentGenomePairs %>% filter(patient_id==patient_id_2)

# Filter to pairs of genomes which are 'almost identical' 
# 177 of our isolates share a match with at least one other isolate from the same subject 
# (by 'almost identical' or >%99.99 identity standards, e.g. fewer than 118 SNPs different)
AlmostIdenticalPairs = Melted_DifferentGenomePairs_SameSubj %>% filter(Almost_Identical)
length(unique(AlmostIdenticalPairs$Genome1))


# A separate one for those which have a 0 SNP distance (identical core genomes)
# 95 of them 
IdenticalPairs = Melted_DifferentGenomePairs_SameSubj %>% filter(Identical)
unpairedisolates_identical = setdiff(unique(Melted$Genome1), unique(IdenticalPairs$Genome1))



# Assign 'representative genome' arbitrarily (even though we'll actually choose the rep based on coverage etc.
# once we have the cluster) One for <5 SNPs (slightly less stringent than requiring 0 distance-- might allow for false positives)
AlmostAlmostIdenticalPairs = Melted_DifferentGenomePairs_SameSubj %>% filter(LessThan5_Snps)
length(unique(AlmostAlmostIdenticalPairs$Genome1))
unpaired = setdiff(unique(Melted$Genome1), unique(IdenticalPairs$Genome1))

# "Clusters" based on "almost almost identical (<5 SNPs)" core genomes 
# AlmostAlmostIdenticalClusters <- list()
# Unique_Genome_List = (unique(AlmostAlmostIdenticalPairs$Genome1))
# Number_AlreadyUnique_Genomes = 221-length(Unique_Genome_List)
# while(length(Unique_Genome_List) > 0){
#   RepGenome = Unique_Genome_List[1]
#   LittleList=(AlmostAlmostIdenticalPairs %>% filter(Genome1==RepGenome))$Genome2
#   AlmostAlmostIdenticalClusters[[RepGenome]] <- LittleList
#   Unique_Genome_List = setdiff(Unique_Genome_List, LittleList)
#   Unique_Genome_List = setdiff(Unique_Genome_List, c(RepGenome))
#   
# }
# # Numbers of isolates we'd have if we subsetted this way 
# Number_AlreadyUnique_Genomes + length(AlmostAlmostIdenticalClusters)

# 'Clusters' based on >99.99% identity core genomes ('Almost Identical')
########################################################################
AlmostIdenticalClusters <- list()

# Unique genome names which are part of pairs
AlmostUnique_Genome_List = (unique(AlmostIdenticalPairs$Genome1))

# Completely unpaired genomes (ones that just need to be sequenced on their own)
Number_AlreadyUnique_Genomes = 221-length(Unique_Genome_List)

unpaired_almostidentical = setdiff(unique(Melted$Genome1), unique(AlmostIdenticalPairs$Genome1))
write.table(unpaired_almostidentical, file="CompletelyUnique_99.99pctlevel.txt", row.names=F, col.names=F, sep="")

# Assign within-subject isolates to clusters 
j=1
length(unique(AlmostIdenticalPairs$patient_id)) # If we just grouped by subject and did no further grouping by identity 
while(length(AlmostUnique_Genome_List) > 0){
  RepGenome = AlmostUnique_Genome_List[1]
  
  LittleList=(AlmostIdenticalPairs %>% filter(Genome1==RepGenome))$Genome2
  AlmostIdenticalClusters[[j]] <- append(LittleList, RepGenome)
  j=j+1
  AlmostUnique_Genome_List = setdiff(AlmostUnique_Genome_List, LittleList)
  AlmostUnique_Genome_List = setdiff(AlmostUnique_Genome_List, c(RepGenome))
}
# Numbers of isolates we'd have if we subsetted this way: 86
Number_AlreadyUnique_Genomes + length(AlmostIdenticalClusters)

# Should maybe add a factor for clusters? 
AlmostIdenticalClusters

# 
##########################
clusterassign <- function(x, cluster_assignment_list) {
    returnval=0
    for(i in (1:(length(cluster_assignment_list)))){
      if(x %in% cluster_assignment_list[[i]]){
        returnval = i
      }
    }
    return(returnval)
}

genomeDF = MappingSubjectTimepoint_order# %>% select(c(Genome1, patient_id, visit,subject_timepoint,visit_color))

genomeDF = genomeDF %>% rowwise() %>%  mutate(cluster = clusterassign(Genome1, AlmostIdenticalClusters))
genomeDF_NonUnique = genomeDF %>% filter(cluster!=0)
# Basically, all subjects just have one cluster except for 124 and 176, which each have 2. 

View((genomeDF_NonUnique %>% group_by(cluster, patient_id)) %>% tally())
# Incorporating phenotypic distributions (identifying
# isolates with equivalent distributions of phenotypes )
# This one is hard: 
# Basically, in my understanding, ARM has taken 3 technical 
# replicates for each of 3 or 9(???) biological replicates for each 
# primary isolate. She averaged the technical replicates, 
# and then took the SD of the biological replicates which is listed here. 

# For a given phenotype, we know:

# The observed values of each isolate's phenotype 
# We know the variability associated with each isolate's phenotype

# One way to decide 'are these equivalent isolates' would be to think about it in terms of:

# Assume the isolates are different & their samples are therefore representative of different true distributions 
# What's the probability that

# The probability youd see each of these samples given that they came from equivalent distributions 

# look at P(their values would be THIS different from one another | they are from the same strain)
# We know the SD associated with each isolate's phenotypic behavior 
#     and P(their values would be THIS different from one another | they are from different strains)

# We know the observed values of each isolate's phenotype,
# and we know the SD of each sample (set of measurements from each isolate)
# 
# In the Kolmogorov-Smirnov test: 
# "Quantifies a distance between the empirical distribution functions  of two samples"
# 

# My original line of thinking was that we should also choose isolate
# 'groups' to choose one rep based on the similarity of their phenotypes.
# 118 SNPs, however, is still a lot to consider those 'the same' at a 
# core genome level(beyond to just choose genomes which are similar that 
# we might be able to use the same long-read scaffold)

# So, the following phenotype analyses are performed only on 'clusters' of *identical* core genomes
####################################################################################################
Normalized_Phenotypes = read.csv("Phenotypes091719.csv")

Normalized_Phenotypes$Genome1 = paste0("DORN",Normalized_Phenotypes$DORN)

PhenotypeMeans = Normalized_Phenotypes %>% select(Genome1, XanthinPhenotype, BiofilmPhenotype)

genomeDF = genomeDF %>% left_join(PhenotypeMeans, by="Genome1")

# "Clusters" based on identical genomes
###################################
IdenticalPairsNonMatches = IdenticalPairs %>% filter(Genome1 != Genome2)
Unique_Genome_List = (unique(IdenticalPairsNonMatches$Genome1))
Number_AlreadyUnique_Genomes = 221-length(Unique_Genome_List)

IdenticalClusters = list()
j=1
length(unique(IdenticalPairsNonMatches$patient_id)) # If we just grouped by subject and did no further grouping by identity 

# Build clusters
#################
while(length(Unique_Genome_List) > 0){
  RepGenome = Unique_Genome_List[1]
  print((IdenticalPairsNonMatches %>% filter(Genome1==RepGenome))$Genome2)
  LittleList=(IdenticalPairsNonMatches %>% filter(Genome1==RepGenome))$Genome2
  IdenticalClusters[[j]] <- append(LittleList, RepGenome)
  j=j+1
  Unique_Genome_List = setdiff(Unique_Genome_List, LittleList)
  Unique_Genome_List = setdiff(Unique_Genome_List, c(RepGenome))
}

# Assign Clusters
##################
genomeDF_uniquepairs = MappingSubjectTimepoint_order
genomeDF_uniquepairs = genomeDF_uniquepairs %>% rowwise() %>%  mutate(cluster = clusterassign(Genome1, IdenticalClusters))
genomeDF_paired = genomeDF_uniquepairs %>% filter(cluster!=0)


UgenomeDF = genomeDF_paired %>% left_join(PhenotypeMeans, by="Genome1")

# Break it up into 4 groups (with 8 clusters displayed in each for ease of visualization)
UgenomeDF1 = UgenomeDF %>% filter(cluster %in% 1:8)
UgenomeDF2 = UgenomeDF %>% filter(cluster %in% 9:16)
UgenomeDF3 = UgenomeDF %>% filter(cluster %in% 17:24)
UgenomeDF4 = UgenomeDF %>% filter(cluster %in% 25:32)

library(grid)
up1 <- ggplot2::ggplot(data=UgenomeDF1, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Clusters 1-8") + guides(col=guide_legend("Cluster"))
up2 <- ggplot2::ggplot(data=UgenomeDF2, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Clusters 9-16") + guides(col=guide_legend("Cluster"))
up3 <- ggplot2::ggplot(data=UgenomeDF3, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Clusters 17-24") + guides(col=guide_legend("Cluster"))
up4 <- ggplot2::ggplot(data=UgenomeDF4, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Clusters 25-32") + guides(col=guide_legend("Cluster"))
postscript("phenotypes.eps")
gridExtra::grid.arrange(up1, up2, up3, up4, top=grid::textGrob("Staphyloxanthin & Biofilm Phenotypes in Clusters of Identical Core Genomes",gp=gpar(fontsize=20,font=3)) )
dev.off()
# Numbers of isolates we'd have if we subsetted this way: 86
Number_AlreadyUnique_Genomes + length(AlmostIdenticalClusters)
#length(Clusters) + Number_AlreadyUnique_Genomes
# 158 isolates after subsetting by zero-SNP core genome distance, shared subject



genomeDF1 = genomeDF %>% filter(cluster %in% 1:5)
genomeDF2 = genomeDF %>% filter(cluster %in% 6:10)
genomeDF3 = genomeDF %>% filter(cluster %in% 11:15)
genomeDF4 = genomeDF %>% filter(cluster %in% 16:20)
newpalette_1 = c("#9900FF", "#FF9900", "#0072B2","#33CC33", "#4D4D4D","#D119A3", "#56B4E9","#CC99FF",
                 "#F0E442", "#CC0000",  "#006600", "#D119A3", "#56B4E9", "#6B24B2", "#339966", "#999999", "#B85C00")
ggplot2::ggplot(data=genomeDF1, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Isolate Clusters 1-5: Biofilm vs. Xanthin Phenotype")
ggplot2::ggplot(data=genomeDF2, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Isolate Clusters 6-10: Biofilm vs. Xanthin Phenotype")
ggplot2::ggplot(data=genomeDF3, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Isolate Clusters 11-15: Biofilm vs. Xanthin Phenotype")
ggplot2::ggplot(data=genomeDF4, aes(x=XanthinPhenotype, y=BiofilmPhenotype, color=as.factor(cluster))) + geom_point() + scale_color_manual(values=newpalette_1) + ggtitle("Isolate Clusters 16-20: Biofilm vs. Xanthin Phenotype")


sd(Normalized_Phenotypes$BiofilmPhenotype)
sd(Normalized_Phenotypes$XanthinPhenotype)
