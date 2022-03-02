
# Amy Campbell
# figuring out what data we have and what we need
# 
library(stringr)
library(stats)
# library(HDMD) for other distance metrics w/ groups
library(tidyr)
library(dplyr)
library(rdist)
library(ggplot2)
library(ggdendro)
library(cluster)
setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/")

# Functions
############
MinMaxNorm = function(x) {
  print(min(x))
  print(max(x))
  return((x- min(x))/(max(x) - min(x)))}


# (1) Identifying which phenotypes I have for which patients
#############################################################
# Mapping from ARM to DORN numbers
ConversionARMS = read.csv("data/AmeliaStrains.csv")
AllPhenotypes=read.csv("data/phenotype_variation_11.06.20.csv")
AllPhenotypes$DORN = paste0("DORN", AllPhenotypes$DORN)
AllPhenotypes = AllPhenotypes %>% dplyr::select(strain, DORN, siderophore, hemolysis, biofilm, xanthin, kinase)


ConversionARMS = ConversionARMS %>% dplyr::select(Strain.Number, doern.lab.bank.., patient, visit)
colnames(ConversionARMS) = c("ARM", "DORN", "Patient", "Visit")

ConversionARMS$DORN= paste0("DORN", ConversionARMS$DORN)

ListSequences = read.csv("mappings/FinalIsolatesList.txt", header=F)


DORNS = data.frame(sapply(ListSequences$V1, function(x) stringr::str_remove(x, "_Final.fasta")))
colnames(DORNS) = c("DORN")


DORNS = DORNS %>% left_join(ConversionARMS, by="DORN")
DORNSPhenotypes = DORNS %>% left_join(AllPhenotypes, by="DORN")
View(DORNSPhenotypes)

set.seed(19104)


justxanthin = read.csv("data/staphyloxanthin_averages_1.13.21.csv")
justxanthin$DORN = paste0("DORN", justxanthin$X)
justxanthin$StaphyloxanthinAverageJan21 = justxanthin$Average
justxanthin = justxanthin %>% dplyr::select(DORN, StaphyloxanthinAverageJan21)
DORNSPhenotypes = DORNSPhenotypes %>% left_join(justxanthin, by="DORN")


phenotypesJan20 = read.csv("data/Phenotypes_01.15.20.csv")
phenotypesJan20$DORN = paste0("DORN", phenotypesJan20$DORN)
phenotypesJan20 = phenotypesJan20 %>% dplyr::select(DORN,XanthinPhenotype, BiofilmPhenotype)
colnames(phenotypesJan20) = c("DORN", "XanthinJan20", "BiofilmJan20")
DORNSPhenotypes = DORNSPhenotypes %>% left_join(phenotypesJan20, by="DORN")


# (2) Clustering by phenotype to identify phenotypically unique subset of isolates
###################################################################################

# Subset to 3 most complete phenotypes
DORNSPhenotypesDistance = DORNSPhenotypes %>% dplyr::select(c("DORN","StaphyloxanthinAverageJan21", "BiofilmJan20", "kinase"))
DORNSPhenotypesDistance = DORNSPhenotypesDistance %>% drop_na()
savenames = DORNSPhenotypesDistance$DORN
DORNSPhenotypesDistance$DORN =NULL

DORNSPhenotypesDistance$StaphyloxanthinAverageJan21 = sapply(DORNSPhenotypesDistance$StaphyloxanthinAverageJan21, function(x) as.numeric(as.character(x)))
DORNSPhenotypesDistance$BiofilmJan20 = sapply(DORNSPhenotypesDistance$BiofilmJan20, function(x) as.numeric(as.character(x)))
DORNSPhenotypesDistance$kinase = sapply(DORNSPhenotypesDistance$kinase, function(x) as.numeric(as.character(x)))
mat= as.matrix(DORNSPhenotypesDistance)

# Min-max (0 to 1) normalize within each phenotype so euclidean distances aren't skewed by scale
DORNSPhenotypesDistanceMat = apply(mat,2,MinMaxNorm )

PhenoTypedistances = pdist(DORNSPhenotypesDistanceMat, metric = "euclidean")

colnames(PhenoTypedistances) = savenames
rownames(PhenoTypedistances) = savenames

# Melt the distance matrix to allow for easy plotting 
Melted = reshape2::melt(PhenoTypedistances)
colnames(Melted) =c("DORN", "DORNCompare", "EuclideanDist")
DORNSPhenotypesPatientMap = DORNSPhenotypes %>% dplyr::select(Patient, DORN)
DORNSPhenotypesPatientMap1 = DORNSPhenotypes %>% dplyr::select(Patient, DORN)
colnames(DORNSPhenotypesPatientMap1) = c("PatientCompare", "DORNCompare")
Melted = Melted %>% left_join(DORNSPhenotypesPatientMap1, by="DORNCompare")
Melted = Melted %>% left_join(DORNSPhenotypesPatientMap, by="DORN")


Melted$PatientDORN = paste(Melted$Patient, Melted$DORN, sep="_")


# Make initial geom_tile() plot and sort by patient
Melted$DORN = factor(Melted$DORN)
Melted$DORNCompare = factor(Melted$DORNCompare)
Melted$comparison = paste(Melted$DORN, Melted$DORNCompare, sep="_")


TilePlot = ggplot(data = Melted, aes(x=DORNCompare, y=DORN, fill=EuclideanDist)) + 
  geom_tile() + scale_fill_gradient(high="white", low="blue") #+ scale_y_discrete(limits=Melted$PatientDORN)

sorted =  data.frame(TilePlot$data) %>% arrange(Patient, DORN)
levelsgood = unique(sorted$DORN)

sorted$Patient
TilePlot$data$DORN = factor(TilePlot$data$DORN, levels=levelsgood)

TilePlot$data$DORNCompare = factor(TilePlot$data$DORNCompare, levels=levelsgood)

dornlevels = data.frame(levelsgood)
colnames(dornlevels) = c("DORN")
dornlevels = dornlevels %>% left_join(DORNSPhenotypesPatientMap, by="DORN")
patientlabels = dornlevels$Patient

colorlist = c("grey14", "grey85")
previous=105
tickcolors=c()
totalindex=1

for (i in patientlabels){
  # 0 if even, 1 if odd
  
  # If still on the same patient # as the previous, don't change anything
  # Otherwise, update previous to the new patient # and increment 
  if(i==previous){
    totalindex = totalindex
  } else{
    totalindex=(totalindex+1)
    previous=i
  }
  index = (totalindex %% 2) + 1
  tickcolors=append(tickcolors,colorlist[index])
}

# save plot with tick marks coded by changes in patient along axis 
newtile = TilePlot + theme(axis.ticks = element_line(colour =tickcolors, size = 4.5), axis.text.x=element_text(angle=270, vjust = .5, hjust=0))
ggsave(newtile, file="PhenotypeEucDistance.pdf", height=30,width=32)

# Hierarchical clustering (average /UPGMA)
distancemat = as.dist(PhenoTypedistances, diag = TRUE)
cresult = hclust(distancemat)
plot(cresult)
orderClust = cresult$labels[cresult$order]
cresultAverage = hclust(distancemat, method="average")
cresultAverageOrder = cresultAverage$labels[cresultAverage$order]
levelsAverage = unique(cresultAverageOrder)
TilePlotAvg = TilePlot
TilePlotAvg$data$DORN = factor(TilePlotAvg$data$DORN, levels=levelsAverage)
TilePlotAvg$data$DORNCompare = factor(TilePlotAvg$data$DORNCompare, levels=levelsAverage)
dendresultAvg = as.dendrogram(cresultAverage)
dataDendro <- dendro_data(dendresultAvg, type = "rectangle")
pdendroAverage = ggdendrogram(dendresultAvg, rotate = FALSE, size = 2)

# Reorder the labels on the heatmap by position on the hierarchical clustering tree
levelsgoodClust = unique(orderClust)
TilePlot$data$DORN = factor(TilePlot$data$DORN, levels=levelsgoodClust)
TilePlot$data$DORNCompare = factor(TilePlot$data$DORNCompare, levels=levelsgoodClust)


dendresult = as.dendrogram(cresult)
dataDendro <- dendro_data(dendresult, type = "rectangle")
pdendro = ggdendrogram(dendresult, rotate = FALSE, size = 2)


SixClusters = cutree(cresult, 6)
dataframedendro = data.frame(cluster = cutree(cresult, 6),
                             states = factor(cresult$labels, levels = cresult$labels[cresult$order]))
pdendro1 = ggdendrogram(dendresult) +
  coord_cartesian(xlim = c(0, nrow(dataframedendro) + 1), ylim = c(0, max(dataDendro$segments$y)), expand = F)  + theme(axis.text.x=element_blank(),axis.text.y=element_blank(),plot.margin=margin( b = -1,t=0,l=0,r=0))  +labs(x=NULL) + geom_hline(yintercept=.75, color="red")

  
TilePlotClust = TilePlot + theme(axis.text.x=element_text(angle=270, vjust = .5, hjust=0), plot.margin=margin(t=-1,b=0,l=0,r=0)) +  coord_cartesian(xlim = c(0, nrow(dataframedendro) + 1), expand = F)
pdendro1

arrangedclusters = ggarrange(pdendro1, TilePlotClust,ncol=1, heights=c(.3,1))

ggsave(arrangedclusters, file="PhenotypeEucDistanceHClust.pdf", height=38,width=32)
ggarrange(pdendro, TilePlotClust, ncol=1)



# K-means clustering by phenotypes
set.seed(19104)

SumSqList = c()
SilScoreList = c()

distmat = dist(DORNSPhenotypesDistanceMat)

for (k in 2:20){
  KMresults = kmeans(DORNSPhenotypesDistanceMat, k, iter.max = 20, nstart = 25)  
  SumSqList = append(SumSqList, KMresults$tot.withinss)
  silscore = cluster::silhouette(KMresults$cluster, dist(DORNSPhenotypesDistanceMat, method="euclidean"))
  SilScoreList = append(SilScoreList, mean(silscore[,3]))
}
gapstatistics = clusGap(DORNSPhenotypesDistanceMat, FUN = kmeans, nstart = 1,
                        K.max = 20, B = 50)

SumSqDF = data.frame(SumSqList)
colnames(SumSqDF) = c("SumSquares")
SumSqDF$k=2:20
SumSqDF$silhouette = SilScoreList

gapDF = data.frame(gapstatistics$Tab)
gapDF$X=row.names(gapDF)
gapDF$k = sapply(gapDF$X, function(x) as.numeric(as.character(x)))

# 4 clusters based on gap stat and elbow
sumsqplot = ggplot(SumSqDF, aes(x=k, y=SumSquares)) + geom_line() + scale_x_discrete(limits=c(2:20)) + ggtitle("Within-Cluster Sum of Squares(Elbow)")
SilhouettePlot = ggplot(SumSqDF, aes(x=k, y=silhouette)) + geom_line() + scale_x_discrete(limits=c(2:20)) + ggtitle("Silhouette Score")
gapplot = ggplot(gapDF, aes(x=k, y=gap)) + geom_line() + scale_x_discrete(limits=2:20) + ggtitle("Gap Statistic")

cowplot::plot_grid(sumsqplot, SilhouettePlot, gapplot, ncol=3)

TilePlotKmeans = TilePlot

KMresults4 = kmeans(DORNSPhenotypesDistanceMat, 4, iter.max = 20, nstart = 25)  

TilePlotKmeans$data$clusterKmeans= KMresults4$cluster

SortedK4 = TilePlotKmeans$data %>% arrange(clusterKmeans)
levelsKmeans4 = unique(SortedK4$DORN)
TilePlotKmeans$data$DORN = factor(TilePlotKmeans$data$DORN, levels=levelsKmeans4)

TilePlotKmeans$data$DORNCompare = factor(TilePlotKmeans$data$DORNCompare, levels=levelsKmeans4)

sort(TilePlotKmeans$data$clusterKmeans)
KmeanColors = TilePlotKmeans$data  %>% filter(DORNCompare=="DORN1000") %>% arrange(clusterKmeans) %>% mutate(colorclust=case_when(clusterKmeans==1 ~ "#FD7446FF", 
                                                                  clusterKmeans==2 ~ "#FED439FF",
                                                                  clusterKmeans==3 ~ "#075149FF", 
                                                                  clusterKmeans==4 ~ "#C80813FF"))

TilePlotKmeans = TilePlotKmeans +  theme(axis.ticks = element_line(colour =KmeanColors$colorclust, size = 4.5), axis.text.x=element_text(angle=270, vjust = .5, hjust=0))  + ggtitle("K-means clustering (K=4)")
ggsave(TilePlotKmeans, file="KMeansEuclideanDistHeatmap.pdf", height=23, width=25)

uniqueRows = (SortedK4 %>% dplyr::select(Patient, clusterKmeans) %>% unique())

BySubj = uniqueRows %>% dplyr::select(Patient, clusterKmeans)
table(BySubj$Patient)


sort(TilePlotKmeans$data$clusterKmeans)
KmeanColors = TilePlotKmeans$data  %>% filter(DORNCompare=="DORN1000") %>% arrange(clusterKmeans) %>% mutate(colorclust=case_when(clusterKmeans==1 ~ "#FD7446FF", 
                                                                                                                                  clusterKmeans==2 ~ "#FED439FF",
                                                                                                                                  clusterKmeans==3 ~ "#075149FF", 
                                                                                                                                  clusterKmeans==4 ~ "#C80813FF"))

TilePlotKmeans = TilePlotKmeans +  theme(axis.ticks = element_line(colour =KmeanColors$colorclust, size = 4.5), axis.text.x=element_text(angle=270, vjust = .5, hjust=0))  + ggtitle("K-means clustering (K=4)")
ggsave(TilePlotKmeans, file="KMeansEuclideanDistHeatmap.pdf", height=23, width=25)



# Just looking what K=3 looks like 
KMresults3 = kmeans(DORNSPhenotypesDistanceMat, 3, iter.max = 20, nstart = 1)  
TilePlotKmeans3 = TilePlotKmeans
TilePlotKmeans3$data$clusterKmeans= KMresults3$cluster

SortedK3 = TilePlotKmeans3$data %>% arrange(clusterKmeans)
levelsKmeans = unique(SortedK3$DORN)
TilePlotKmeans3$data$DORN = factor(TilePlotKmeans3$data$DORN, levels=levelsKmeans)
TilePlotKmeans3$data$DORNCompare = factor(TilePlotKmeans3$data$DORNCompare, levels=levelsKmeans)
TilePlotKmeans3 = TilePlotKmeans3 + ggtitle("K-means clustering (K=3)")
sort(TilePlotKmeans$data$clusterKmeans)
uniqueRows = (SortedK5 %>% dplyr::dplyr::select(Patient, clusterKmeans) %>% unique())

uniqueRows3 = (SortedK3 %>% dplyr::select(Patient, clusterKmeans) %>% unique())
KmeanColors3 = TilePlotKmeans3$data  %>% filter(DORNCompare=="DORN1000") %>% arrange(clusterKmeans) %>% mutate(colorclust=case_when(clusterKmeans==1 ~ "#FD7446FF", 
                                                                                                                                  clusterKmeans==2 ~ "#FED439FF",
                                                                                                                                  clusterKmeans==3 ~ "#075149FF", 
                                                                                                                                  ))
TilePlotKmeans3 = TilePlotKmeans3 + theme(axis.ticks = element_line(colour =KmeanColors3$colorclust, size = 4.5), axis.text.x=element_text(angle=270, vjust = .5, hjust=0))  + ggtitle("K-means clustering (K=3)")
ggsave(TilePlotKmeans3, file="KMeansEuclideanDistHeatmap_K3.pdf", height=23, width=25)

