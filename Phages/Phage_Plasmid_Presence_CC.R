library("ggstar")
library("ggplot2")
library("cowplot")
library("stringr")
library(dplyr)
genes = read.csv("Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence.csv")
phage_plasmid_pres = read.csv("Documents/DataInputGithub/data/IntraPatient/GenePresence_ByPhage_Plasmid_Updated_Hclust.csv")
CCs = read.csv("Documents/DataInputGithub/data/Phylogeny2022Data/CCMapPlotting.csv")
CCs$Genome = CCs$DORN
CCs = CCs %>% select(Genome, CCLabel)
phage_presence_absence = read.csv("Documents/DataInputGithub/data/IntraPatient/Phages/PresenceAbsenceHClustPhages.csv")
plasmid_presence_absence = read.csv("Documents/DataInputGithub/data/IntraPatient/Plasmids/Filtered_Plasmid_Presence_Absence.csv")
IsolateInfo = read.csv("Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")
IsolateInfo$Genome = paste0("DORN", IsolateInfo$Doern.lab.bank.)

# EAP wansn't in the representative sequence for any phages or plasmids 
########################################################################
EAP = phage_plasmid_pres %>% filter(Annotation=="extracellular adherence protein Eap/Map")
rowSums(EAP[3:(ncol(EAP) -2)])

# Neither was deoC
###################
deoC = phage_plasmid_pres %>% filter(grepl(x=Annotation, "deoxyribose-phosphate"))
rowSums(deoC[3:(ncol(deoC) -2)])

######################
# Figure S1: Plasmids
######################
plasmid_presence_absence  = plasmid_presence_absence %>% select(X, colnames(plasmid_presence_absence)[colnames(plasmid_presence_absence) %in% colnames(phage_plasmid_pres)])


plasmid_presence_absence$Genome = plasmid_presence_absence$X
plasmid_presence_absence$X.1=NULL
plasmid_presence_absence$X = NULL
plasmid_presence_absence = plasmid_presence_absence %>% left_join(IsolateInfo %>% select(Genome, patient_id), by="Genome")
plasmid_presence_absence = plasmid_presence_absence %>% left_join(CCs,by="Genome")
PatientsPerPlasmidCC=data.frame()


for(plasmidNum in 1:(ncol(plasmid_presence_absence)-3)){
  

  plasmidDF = plasmid_presence_absence[,c(plasmidNum, ncol(plasmid_presence_absence), ncol(plasmid_presence_absence)-1)]
  
  colnames(plasmidDF) = c("Plasmid","CC","Patient")

  NumDFUs_With_Plasmid_CC=data.frame(table(plasmidDF %>% filter(Plasmid==1) %>% unique() %>% select(CC)))
  colnames(NumDFUs_With_Plasmid_CC) = c("CC", "NumPatients")
  AllCCs = CCs %>% select(CCLabel) %>% unique()
  colnames(AllCCs) = "CC"
  NumDFUs_With_Plasmid_CC = AllCCs %>% left_join(NumDFUs_With_Plasmid_CC, by="CC")
  NumDFUs_With_Plasmid_CC$Plasmid = colnames(plasmid_presence_absence)[plasmidNum]
  NumDFUs_With_Plasmid_CC[is.na(NumDFUs_With_Plasmid_CC)] <- 0
  
  
  PatientsPerPlasmidCC = rbind(PatientsPerPlasmidCC, NumDFUs_With_Plasmid_CC)

}
colnames(PatientsPerPlasmidCC) = c("CC", "NumPatients", "Plasmid")

sort(unique(PatientsPerPlasmidCC$NumPatients))


PatientsPerPlasmidCC = PatientsPerPlasmidCC %>% mutate(NumDFUcategory = case_when(NumPatients==0 ~ "0", 
                                                                              NumPatients==1 ~ "1", 
                                                                              NumPatients %in% 2:7 ~ "2",
                                                                              NumPatients %in% 8:14 ~ "3"))

# Cluster plasmids by their gene content 
plasmid_Gene_pres = phage_plasmid_pres %>% select(colnames(phage_plasmid_pres)[!grepl("HclustPhage", colnames(phage_plasmid_pres))])
savegenes = plasmid_Gene_pres$Gene
plasmid_Gene_pres$Gene=NULL
plasmid_Gene_pres$X.1=NULL
plasmid_Gene_pres$X=NULL
plasmid_Gene_pres$Annotation=NULL

SavePlasmidNames = colnames(plasmid_Gene_pres)
GenePresAbsDist = dist(t(plasmid_Gene_pres), method="binary")

clusteringPlasmidsByGenes = hclust(GenePresAbsDist, method="average")

OrderPlasmidsByGeneContent = SavePlasmidNames[clusteringPlasmidsByGenes$order]
plasmid_presence_absence$patient_id=NULL

plasmid_CC = plasmid_presence_absence
plasmid_CC$Genome = NULL

plasmid_CC = plasmid_CC %>% group_by(CCLabel) %>% summarize_all(max)


plasmid_CC_melt = plasmid_CC %>% reshape2::melt(id.vars=c("CCLabel"))
colnames(plasmid_CC_melt) = c("CC", "Plasmid", "PresAbs")

plasmid_CC_melt = plasmid_CC_melt %>% left_join(PatientsPerPlasmidCC,by=c("Plasmid", "CC"))


ccorder = c("CC8", "CC15", "CC20", "CC1","CC97","CC7","CC12","CC72","CC9","CC5",
            "CC59","CC45","CC30")


plasmidvals = c("white","skyblue1", "mediumblue","navyblue")
# shapes =  starshapeobj=c(11,15)
plasmid_plot = ggplot(plasmid_CC_melt, aes(x=CC, y=Plasmid, fill=factor(NumDFUcategory), color=factor(NumDFUcategory))) + geom_star(starshape=15, size=2) +
  scale_color_manual(values=plasmidvals) + scale_fill_manual(values=plasmidvals)+ 
  theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), 
        axis.text.y=element_text(size=2),
        legend.position="none", axis.line.x=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        panel.background = element_blank(), axis.title.y=element_blank(),
        axis.title.x=element_blank())+scale_x_discrete(position = "top") 


plasmid_plot$data$Plasmid = factor(plasmid_plot$data$Plasmid, levels=OrderPlasmidsByGeneContent)
plasmid_plot$data$CC = factor(plasmid_plot$data$CC, levels=ccorder)

# plasmid_Gene_pres$rowsums = rowSums(plasmid_Gene_pres)
# plasmid_Gene_pres$Gene = savegenesplasmid_Gene_pres
# plasmid_Gene_pres = plasmid_Gene_pres %>% filter(rowsums>0)
# savegenessubset=plasmid_Gene_pres$Gene
# plasmid_Gene_pres$rowsums=NULL
# plasmid_Gene_pres$Gene = NULL
# 
# genedists = dist(plasmid_Gene_pres, method="binary")
# 
# savegenessubset = savegenessubset[(hclust(genedists, method="average"))$order]
# 
# plasmid_Gene_pres$Gene = savegenessubset
# 
# MeltedGenesPlasmid = reshape2::melt(plasmid_Gene_pres,id.vars=c("Gene"))

# Plot phages by CC  
####################
phage_presence_absence$X=NULL
GenomeSave=phage_presence_absence$Genome
phage_presence_absence = phage_presence_absence %>% left_join(IsolateInfo %>% select(Genome, patient_id), by="Genome")
phage_presence_absence = phage_presence_absence %>% left_join(CCs,by="Genome")


phage_presence_absence$Genome=NULL
PatientsPerPhageCC=data.frame()

for(phageNum in 1:(ncol(phage_presence_absence)-2)){
  
  
  phageDF = phage_presence_absence[,c(phageNum, ncol(phage_presence_absence), ncol(phage_presence_absence)-1)]
  
  colnames(phageDF) = c("Phage","CC","Patient")
  
  NumDFUs_With_phage_CC=data.frame(table(phageDF %>% filter(Phage==1) %>% unique() %>% select(CC)))
  colnames(NumDFUs_With_phage_CC) = c("CC", "NumPatients")
  AllCCs = CCs %>% select(CCLabel) %>% unique()
  colnames(AllCCs) = "CC"
  NumDFUs_With_phage_CC = AllCCs %>% left_join(NumDFUs_With_phage_CC, by="CC")
  NumDFUs_With_phage_CC$Phage = colnames(phage_presence_absence)[phageNum]
  NumDFUs_With_phage_CC[is.na(NumDFUs_With_phage_CC)] <- 0
  
  
  PatientsPerPhageCC = rbind(PatientsPerPhageCC, NumDFUs_With_phage_CC)
  
}
colnames(PatientsPerPhageCC) = c("CC", "NumPatients", "Phage")

max(PatientsPerPhageCC$NumPatients)

PatientsPerPhageCC = PatientsPerPhageCC %>% mutate(NumDFUcategory = case_when(NumPatients==0 ~ "0", 
                                                                                  NumPatients==1 ~ "1", 
                                                                                  NumPatients %in% 2:7 ~ "2",
                                                                                  NumPatients %in% 8:13 ~ "3"))


phage_Gene_pres = phage_plasmid_pres %>% select("Gene", colnames(phage_plasmid_pres)[grepl("HclustPhage", colnames(phage_plasmid_pres))])
savegenesPhage = phage_Gene_pres$Gene
phage_Gene_pres$Gene=NULL
phage_Gene_pres$X.1=NULL
phage_Gene_pres$X=NULL
phage_Gene_pres$Annotation=NULL

SavePhageNames = colnames(phage_Gene_pres)
GenePresAbsDist = dist(t(phage_Gene_pres), method="binary")

clusteringPhagebyGenes = hclust(GenePresAbsDist, method="average")

OrderPhagesByGeneContent = SavePhageNames[clusteringPhagebyGenes$order]
phage_presence_absence$patient_id=NULL

phage_CC = phage_presence_absence

phage_CC = phage_CC %>% group_by(CCLabel) %>% summarize_all(max)


phage_CC_melt = phage_CC %>% reshape2::melt(id.vars=c("CCLabel"))
colnames(phage_CC_melt) = c("CC", "Phage", "PresAbs")

phage_CC_melt = phage_CC_melt %>% left_join(PatientsPerPhageCC,by=c("Phage", "CC"))


ccorder = c("CC8", "CC15", "CC20", "CC1","CC97","CC7","CC12","CC72","CC9","CC5",
            "CC59","CC45","CC30")

phagevals = c("white","#AFE1AF","#4CBB17","#024B30")

phage_plot = ggplot(phage_CC_melt, aes(x=CC, y=Phage, fill=factor(NumDFUcategory), color=factor(NumDFUcategory))) + geom_star(starshape=11, size=2) +
  scale_color_manual(values=phagevals) + scale_fill_manual(values=phagevals)+ 
  theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), 
        axis.text.y=element_text(size=2),
        legend.position="none", axis.line.x=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        panel.background = element_blank(), axis.title.y=element_blank(),
        axis.title.x=element_blank())+scale_x_discrete(position = "top")  #+ coord_fixed(ratio=1)

phage_plot$data$Phage = factor(phage_plot$data$Phage, levels=OrderPhagesByGeneContent)
phage_plot$data$CC = factor(phage_plot$data$CC, levels=ccorder)




Cluster_All_Genes_PresenceAbsence = phage_plasmid_pres
Cluster_All_Genes_PresenceAbsence = Cluster_All_Genes_PresenceAbsence %>% select(-Annotation,-X, -X.1)

Cluster_All_Genes_PresenceAbsence$sums = rowSums(Cluster_All_Genes_PresenceAbsence[,1:ncol(Cluster_All_Genes_PresenceAbsence)-1])

Cluster_All_Genes_PresenceAbsence = Cluster_All_Genes_PresenceAbsence %>% filter(sums>0)

saveGenesAll = Cluster_All_Genes_PresenceAbsence$Gene
Cluster_All_Genes_PresenceAbsence$Gene = NULL

geneorder = (hclust(dist(Cluster_All_Genes_PresenceAbsence, method="binary")))$order

geneord_string = saveGenesAll[geneorder]
Cluster_All_Genes_PresenceAbsence$Gene = saveGenesAll


# Genes by Plasmid
PlasmidGenes = Cluster_All_Genes_PresenceAbsence %>% select(Gene, colnames(Cluster_All_Genes_PresenceAbsence)[colnames(Cluster_All_Genes_PresenceAbsence)%in% colnames(plasmid_CC) ])

PlasmidGenesMelt = reshape2::melt(PlasmidGenes, id.vars=c("Gene"))
PlasmidGenePlot = ggplot(PlasmidGenesMelt, aes(x=Gene,y=variable,fill=factor(value))) + geom_tile() + scale_fill_manual(values=c("white","black"))
PlasmidGenePlot$data$variable = factor(PlasmidGenePlot$data$variable, levels=OrderPlasmidsByGeneContent)
PlasmidGenePlot$data$Gene = factor(PlasmidGenePlot$data$Gene, levels=geneord_string)


# Genes by phage

PhageGenes = Cluster_All_Genes_PresenceAbsence %>% select(Gene, colnames(Cluster_All_Genes_PresenceAbsence)[colnames(Cluster_All_Genes_PresenceAbsence)%in% colnames(phage_CC) ])

PhageGenesMelt = reshape2::melt(PhageGenes, id.vars=c("Gene"))
PhageGenePlot = ggplot(PhageGenesMelt, aes(x=Gene,y=variable,fill=factor(value))) + geom_tile() + scale_fill_manual(values=c("white","black"))
PhageGenePlot$data$variable = factor(PhageGenePlot$data$variable, levels=OrderPhagesByGeneContent)
PhageGenePlot$data$Gene = factor(PhageGenePlot$data$Gene, levels=geneord_string)
PhageGenePlot = PhageGenePlot + theme(axis.text.x=element_text(angle=90))
ggsave(PhageGenePlot, file="Documents/Saureus_Genomics_Paper/PhageGenesPlot.pdf",width=45,height=5)
ggsave(PlasmidGenePlot + theme(axis.text.x=element_text(angle=90)), file="Documents/Saureus_Genomics_Paper/PlasmidGenesPlot.pdf",width=45,height=5)

PlasmidGenePlot= PlasmidGenePlot + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) #+ theme_classic()
PhageGenePlot= PhageGenePlot + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) #+ theme_classic()

cowplotGenes = cowplot::plot_grid(PlasmidGenePlot + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), PhageGenePlot+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), align = "v",ncol=1, rel_heights = c(39,14))
ggsave(cowplotGenes, file="Documents/Saureus_Genomics_Paper/GenePresenceHeatmapHGTs.pdf",width=10,height=20)
ggsave(plasmid_plot, width=4, height=6, file="Documents/Saureus_Genomics_Paper/justlabels.pdf")

multiplo_obj = cowplot::plot_grid(plasmid_plot + theme( axis.text.x=element_blank(), axis.title.x=element_blank()),phage_plot+theme( axis.text.x=element_blank(), axis.title.x=element_blank()),align = "v",ncol=1,axis="t", rel_heights = c(39,14))
ggsave(multiplo_obj, file="Documents/Saureus_Genomics_Paper/PhagePlasmidPresence_CC.pdf", height=6,width=1.8)



Megaplot = cowplot::plot_grid(multiplo_obj,cowplotGenes,align="h", axis='l', ncol=2, rel_widths=c(1,2))
ggsave(Megaplot, file="Documents/Saureus_Genomics_Paper/PhagePlasmidPresence_genes_CC.pdf", height=6,width=5.5)


# # Order phages by Hclust
# phage_presence_absence$DORN = phage_presence_absence$X
# phage_presence_absence$X.1=NULL
# phage_presence_absence$X=NULL
# phage_CC = phage_presence_absence %>% left_join(CCs,by="DORN")
# phage_presence_absence$DORN=NULL
# phage_CC = phage_CC %>% group_by(CCLabel) %>% summarize_all(max)
# saveCCs = (phage_CC$CCLabel)
# phage_CC$CCLabel = NULL
# phage_CC$DORN = NULL
# TransposedCC_Phages = data.frame(t(phage_CC))
# colnames(TransposedCC_Phages) = saveCCs
# hclustPhages = hclust(dist(TransposedCC_Phages, method="binary"))
# phage_order = row.names(TransposedCC_Phages)[hclustPhages$order]
# phage_CC$CCLabel = saveCCs
# phage_CC_melt = phage_CC %>% reshape2::melt(id.vars="CCLabel")
# 
# phagevals = c("white","#AFE1AF","#4CBB17","#024B30")
# 
# phage_plot = ggplot(phage_CC_melt, aes(x=CCLabel, y=variable, fill=factor(value), color=factor(value))) + geom_star(starshape=11, size=2) +
#   scale_color_manual(values=phagevals) + scale_fill_manual(values=phagevals)+ 
#   theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), 
#         axis.text.y=element_text(size=2),
#         legend.position="none", axis.line.x=element_blank(),
#         axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
#         panel.background = element_blank(), axis.title.y=element_blank(),
#         axis.title.x=element_blank())+scale_x_discrete(position = "top")  + coord_fixed(ratio=1)
# 
# phage_plot$data$variable = factor(phage_plot$data$variable, levels=phage_order)
# phage_plot$data$CCLabel = factor(phage_plot$data$CCLabel, levels=ccorder)
# 
# 
# 
# Phage_plasmid_CC = Phage_plasmid_presence %>% left_join(CCs,by="Genome")
# 
# 
# Phage_plasmid_CC$DORN=NULL
# Phage_plasmid_CC = Phage_plasmid_CC %>% group_by(CCLabel) %>% summarize_all(max)
# saveCCs = (Phage_plasmid_CC$CCLabel)
# Phage_plasmid_CC$CCLabel = NULL
# TransposedCC_hgts = data.frame(t(Phage_plasmid_CC))
# colnames(TransposedCC_hgts) = saveCCs
# 
# 
# 
# multiplo_obj = cowplot::plot_grid(phage_plot, plasmid_plot,align = "v",ncol=1, rel_heights = c(7,5))
# ggsave(multiplo_obj, file="Documents/Saureus_Genomics_Paper/PhagePlasmidPresence_CC.pdf", height=14.1,width=1.8)
# 
# phagelistmashpath = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/mashresults_phages/"
# phagelist_mash = list.files(phagelistmashpath)
# phage_mash_DF = data.frame()
# for(l in phagelist_mash){
#   mashresult = read.csv2(paste0(phagelistmashpath, l), sep="\t")
#   colnames(mashresult) = c("NCBI", "path", "dist", "p", "matching")
#   littledb = mashresult[1,c("NCBI", "matching")]
#   phagename=str_replace(l, ".tsv", "")
#   littledb$Phage=phagename
#   phage_mash_DF = rbind(phage_mash_DF,littledb )
# }
# 
# 
# mashresult = read.csv2("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/mashresults_phages/Phage11.tsv", sep="\t", header=F)
# 
