# Amy Campbell
# Seeing if the 80% phage distance thing makes sense to use 
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(reshape2)
MashTSVpath = "~/Documents/DataInputGithub/data/IntraPatient/Phages/cdhit_80/MashDistsAll1000.tsv"
ANI_identityTSVpath = "~/Documents/DataInputGithub/data/IntraPatient/Phages/cdhit_80/ANIm_percentage_identity.tab"
AllSaureusIncluded = read.csv("~/Documents/DataInputGithub/data/Phylogeny2022Data/CCMapPlotting.csv")
Summary_sequences = read.csv("~/Documents/DataInputGithub/data/IntraPatient/Phages/SummaryViralSequences.csv")
RenameFastasScript = "~/Desktop/GriceLabGit/DFUStrainsWGS/Phages/ReAssignPhageClusters.sh"

Summary_sequences$CommonID = sapply(Summary_sequences$SequenceName, function(x) str_split(x, pattern=":")[[1]][1])
Summary_sequences$FastaIDWithoutAssignment = sapply(Summary_sequences$CommonID, function(x) str_replace(x, "\\|", "_"))

mashtable = read.csv2(MashTSVpath, sep="\t")
ANItable=read.csv2(ANI_identityTSVpath, sep="\t")
row.names(ANItable) = ANItable$X

mashtable$phageinstance1 = sapply(mashtable$phageinstance1,function(x) str_remove_all(x, "\\/home/acampbe\\/DFU\\/data\\/WGS_2020\\/PhageResults\\/CheckVResults_Parsed\\/cdhit_80\\/IndividualFastas\\/"))
mashtable$phageinstance2 = sapply(mashtable$phageinstance2,function(x) str_remove_all(x, "\\/home/acampbe\\/DFU\\/data\\/WGS_2020\\/PhageResults\\/CheckVResults_Parsed\\/cdhit_80\\/IndividualFastas\\/"))
mashtable$phageinstance1 = sapply(mashtable$phageinstance1,function(x) str_remove_all(x, ".fasta"))
mashtable$phageinstance2 = sapply(mashtable$phageinstance2,function(x) str_remove_all(x, ".fasta"))


mashtable$distance = sapply(mashtable$dist, function(x) as.numeric(as.character(x)))
dist_plot = ggplot(mashtable, aes(x=factor(phageinstance1), y=factor(phageinstance2), fill=distance)) + geom_tile() + scale_fill_viridis(option="plasma") 
dist_plot = dist_plot + theme(axis.text.x=element_text(angle=90))
plotorder = unique(sort(dist_plot$data$phageinstance1))

dist_plot$data$phageinstance1 = factor(dist_plot$data$phageinstance1, levels=plotorder)
dist_plot$data$phageinstance2 = factor(dist_plot$data$phageinstance2, levels=plotorder)
ggsave(dist_plot, file="Documents/Saureus_Genomics_Paper/PhageClusters_k21.pdf", width=40,height=40)


mashmat = dcast(data = mashtable,formula = phageinstance1~phageinstance2,value.var = "distance")
hclust(mashmat)
savephageinst = mashmat$phageinstance1
mashmat$phageinstance1 = NULL
hclustresult = hclust(as.dist(mashmat), method="complete")
hclustresult_single = hclust(as.dist(mashmat), method="single")

clusters = (cutree(hclustresult,h=0.06))
clusters_single = (cutree(hclustresult_single,h=0.05))
clusters_pt05 =(cutree(hclustresult,h=0.05))
clusterassign = data.frame(clusters)
orderplot=names(sort(clusters))

dist_plot$data 
dist_plot$data$phageinstance1 = factor(dist_plot$data$phageinstance1, levels=orderplot)
dist_plot$data$phageinstance2 = factor(dist_plot$data$phageinstance2, levels=orderplot)


ggsave(dist_plot, file="Documents/Saureus_Genomics_Paper/PhageClusters_k21_HC.pdf", width=40,height=40)
plot(x = hclustresult, labels =  row.names(hclustresult), cex = 0.5) + abline(.06,0,col="red")

clustersassigned = data.frame(clusters)
clustersassigned$clusterX = clustersassigned$clusters
clustersassigned$clusterY = clustersassigned$clusters

clustersassigned$phageinstance1 = row.names(clustersassigned)
clustersassigned$phageinstance2 = row.names(clustersassigned)

# Plot with new hierarchical cluster labels
###########################################
mashtable = mashtable %>% left_join(clustersassigned %>% select(clusterX,phageinstance1 ), by="phageinstance1")
mashtable = mashtable %>% left_join(clustersassigned %>% select(clusterY,phageinstance2 ), by="phageinstance2")
mashtable$c1 = paste(mashtable$clusterX, mashtable$phageinstance1, sep="_")
mashtable$c2 = paste(mashtable$clusterY, mashtable$phageinstance2, sep="_")

dist_plot_hclust = ggplot(mashtable, aes(x=factor(c1), y=factor(c2), fill=distance)) + geom_tile() + scale_fill_viridis(option="plasma") 
dist_plot_hclust = dist_plot_hclust + theme(axis.text.x=element_text(angle=90))
plotorder = unique( (mashtable %>% arrange(clusterX, c1))$c1)
dist_plot_hclust$data$c1 = factor(dist_plot_hclust$data$c1, levels=plotorder)
dist_plot_hclust$data$c2 = factor(dist_plot_hclust$data$c2, levels=plotorder)
ggsave(dist_plot_hclust, file="Documents/Saureus_Genomics_Paper/PhageClusters_k21_LabeledHClust.pdf", width=40,height=40)



PctIdentity = reshape2::melt(ANItable, id.vars=c('X'))
PctIdentity$value = sapply(PctIdentity$value, function(x) as.numeric(as.character(x)))
orderplot_pct = sapply(orderplot, function(x) str_replace(x, "\\|", "_"))

PctIdentity$OneMinusPctID =1-PctIdentity$value 
PctIdentity_plot = ggplot(PctIdentity, aes(x=factor(X), y=factor(variable), fill=OneMinusPctID)) + geom_tile() + scale_fill_viridis(option="plasma") 
PctIdentity_plot$data$X = factor(PctIdentity_plot$data$X, levels=orderplot_pct)
PctIdentity_plot$data$variable = factor(PctIdentity_plot$data$variable, levels=orderplot_pct)
ggsave(PctIdentity_plot, file="Documents/Saureus_Genomics_Paper/PctID_HclustClusters.pdf", width=40,height=40)


orderplot_cdhit = sapply(plotorder, function(x) str_replace(x, "\\|", "_"))
PctIdentity_plot$data$X = factor(PctIdentity_plot$data$X, levels=orderplot_cdhit)
PctIdentity_plot$data$variable = factor(PctIdentity_plot$data$variable, levels=orderplot_cdhit)
ggsave(PctIdentity_plot, file="Documents/Saureus_Genomics_Paper/PctID_CDhitclusters.pdf", width=40,height=40)


# It seems like actually the hclust with distance cutoff .06  makes the most sense
clusterassign$ID = row.names(clusterassign)
clusterassign$Genome = sapply(clusterassign$ID,function(x) str_split(x, "_")[[1]][2])


colnames(clusterassign) = c("HierarchicalCluster", "FastaID", "Genome")
clusterassign$FastaName = sapply(clusterassign$FastaID, function(x) str_replace(x, "\\|", "_"))
clusterassign$FastaIDWithoutAssignment=sapply(clusterassign$FastaName, function(x) paste(str_split(x,"_")[[1]][2:4], collapse="_"))
clusterassign$FastaName = paste0(clusterassign$FastaName, ".fasta")
clusterassign$NewFastaName = paste0("HclustPhage",clusterassign$HierarchicalCluster,"_", clusterassign$FastaIDWithoutAssignment, ".fasta")
clusterassign$OldFastaPath = paste0("/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/IndividualFastas/", clusterassign$FastaName)
clusterassign$NewFastaPath = paste0("/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/IndividualFastas_Hclust/", clusterassign$NewFastaName)


clusterassign_presenceAbsence = clusterassign %>% select(HierarchicalCluster, Genome)
clusterassign_presenceAbsence$PhageName = paste0("HclustPhage", clusterassign_presenceAbsence$HierarchicalCluster)
clusterassign_presenceAbsence$Present = 1

AllSaureusIncluded$Genome = AllSaureusIncluded$DORN
clusterassign_presenceAbsence = clusterassign_presenceAbsence %>% select(Genome, PhageName, Present)
rownames(clusterassign_presenceAbsence) = 1:nrow(clusterassign_presenceAbsence)

Presence_Absence_Phages = dcast(data=clusterassign_presenceAbsence,formula=Genome~PhageName,value.var="Present")
Presence_Absence_Phages = AllSaureusIncluded %>% select(Genome) %>% left_join(Presence_Absence_Phages,by="Genome")
Presence_Absence_Phages[is.na(Presence_Absence_Phages)]<-0
write.csv(Presence_Absence_Phages, "~/Documents/DataInputGithub/data/IntraPatient/Phages/PresenceAbsenceHClustPhages.csv")


Lengths = clusterassign %>% select(FastaIDWithoutAssignment, HierarchicalCluster) %>% left_join(Summary_sequences %>% select(FastaIDWithoutAssignment,Length ), by="FastaIDWithoutAssignment")
ggplot(Lengths, aes(x=factor(HierarchicalCluster), y=Length)) + geom_boxplot() + xlab("Hierarchical Phage Cluster") + ylab("Lengths of Sequences Assigned")
Lengths %>% filter(HierarchicalCluster==10) %>% arrange(Length)


RenameScriptCommandlist=c("#!/bin/bash","#Amy Campbell","#Renaming the CD-hit-assigned phages to their hierarchical clustering cluster and making a copy of them that reflects this",
                          "","# Note: this script doesn't alone change the fasta header, so the fasta headers on the new files are still old phage cluster assignments")

for(i in 1:nrow(clusterassign)){
  oldfasta= clusterassign[i,"OldFastaPath"]
  newfasta=clusterassign[i,"NewFastaPath"]
  RenameScriptCommandlist = append(RenameScriptCommandlist, paste0("cp ",oldfasta, " ", newfasta ))
}
writeLines(RenameScriptCommandlist, RenameFastasScript)
