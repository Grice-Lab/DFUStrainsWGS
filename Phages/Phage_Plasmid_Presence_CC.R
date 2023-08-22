library("ggstar")
library("ggplot2")
library("cowplot")


genes = read.csv("Documents/DFUData/RoaryOutputAll/RoaryResultsPGAP2022/gene_presence_absence.csv")
phage_plasmid_pres = read.csv("Downloads/train_ride/GenePresence_ByPhage_Plasmid_Updated.csv")
CCs = read.csv("Documents/DFUData/CCMapPlotting.csv")
# EAP wansn't in the representative sequence for any phages or plasmids 
########################################################################
EAP = phage_plasmid_pres %>% filter(Annotation=="extracellular adherence protein Eap/Map")
rowSums(EAP[3:(ncol(EAP) -2)])

# Neither was deoC
###################
deoC = phage_plasmid_pres %>% filter(grepl(x=Annotation, "deoxyribose-phosphate"))
rowSums(deoC[3:(ncol(deoC) -2)])

phage_presence_absence = read.csv("Downloads/train_ride/PhagePresenceAbsence.csv")
phage_presence_absence = (phage_presence_absence %>% select(colnames(phage_presence_absence)[colnames(phage_presence_absence) %in% colnames(phage_plasmid_pres) ]))

plasmid_presence_absence = read.csv("Downloads/train_ride/Filtered_Plasmid_Presence_Absence.csv")
plasmid_presence_absence = (plasmid_presence_absence %>% select(colnames(plasmid_presence_absence)[colnames(plasmid_presence_absence) %in% colnames(phage_plasmid_pres) ]))

plasmid_presence_absence  = plasmid_presence_absence %>% arrange(X)
phage_presence_absence = phage_presence_absence %>% arrange(X)


Phage_plasmid_presence = phage_presence_absence %>% left_join(plasmid_presence_absence, by="X")
Phage_plasmid_presence$DORN = Phage_plasmid_presence$X
Phage_plasmid_presence$X=NULL
Phage_plasmid_presence$X.1 = NULL

# Order plasmids by Hclust
plasmid_presence_absence$DORN = plasmid_presence_absence$X
plasmid_presence_absence$DORN = plasmid_presence_absence$X
plasmid_presence_absence$X.1=NULL
plasmid_presence_absence$X=NULL
plasmid_CC = plasmid_presence_absence %>% left_join(CCs,by="DORN")
plasmid_CC$DORN=NULL
plasmid_CC = plasmid_CC %>% group_by(CCLabel) %>% summarize_all(max)
saveCCs = (plasmid_CC$CCLabel)
plasmid_CC$CCLabel = NULL
TransposedCC_Plasmids = data.frame(t(plasmid_CC))
colnames(TransposedCC_Plasmids) = saveCCs
hclustPlasmids = hclust(dist(TransposedCC_Plasmids, method="binary"))
plasmid_order = row.names(TransposedCC_Plasmids)[hclustPlasmids$order]
plasmid_CC$CCLabel = saveCCs
plasmid_CC_melt = plasmid_CC %>% reshape2::melt(id.vars="CCLabel")

ccorder = c("CC8", "CC15", "CC20", "CC1","CC97","CC7","CC12","CC72","CC9","CC5",
            "CC59","CC45","CC30")
plasmidvals = c("white", "darkblue")
# shapes =  starshapeobj=c(11,15)
plasmid_plot = ggplot(plasmid_CC_melt, aes(x=CCLabel, y=variable, fill=factor(value), color=factor(value))) + geom_star(starshape=15, size=2) +
  scale_color_manual(values=plasmidvals) + scale_fill_manual(values=plasmidvals)+ 
  theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), 
        axis.text.y=element_text(size=2),
        legend.position="none", axis.line.x=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        panel.background = element_blank(), axis.title.y=element_blank(),
        axis.title.x=element_blank())+scale_x_discrete(position = "top") 


plasmid_plot$data$variable = factor(plasmid_plot$data$variable, levels=plasmid_order)
plasmid_plot$data$CCLabel = factor(plasmid_plot$data$CCLabel, levels=ccorder)


# Order phages by Hclust
phage_presence_absence$DORN = phage_presence_absence$X
phage_presence_absence$X.1=NULL
phage_presence_absence$X=NULL
phage_CC = phage_presence_absence %>% left_join(CCs,by="DORN")
phage_presence_absence$DORN=NULL
phage_CC = phage_CC %>% group_by(CCLabel) %>% summarize_all(max)
saveCCs = (phage_CC$CCLabel)
phage_CC$CCLabel = NULL
phage_CC$DORN = NULL
TransposedCC_Phages = data.frame(t(phage_CC))
colnames(TransposedCC_Phages) = saveCCs
hclustPhages = hclust(dist(TransposedCC_Phages, method="binary"))
phage_order = row.names(TransposedCC_Phages)[hclustPhages$order]
phage_CC$CCLabel = saveCCs
phage_CC_melt = phage_CC %>% reshape2::melt(id.vars="CCLabel")

phagevals = c("white","#024B30")

phage_plot = ggplot(phage_CC_melt, aes(x=CCLabel, y=variable, fill=factor(value), color=factor(value))) + geom_star(starshape=11, size=2) +
  scale_color_manual(values=phagevals) + scale_fill_manual(values=phagevals)+ 
  theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), 
        axis.text.y=element_text(size=2),
        legend.position="none", axis.line.x=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        panel.background = element_blank(), axis.title.y=element_blank(),
        axis.title.x=element_blank())+scale_x_discrete(position = "top")  + coord_fixed(ratio=1)

phage_plot$data$variable = factor(phage_plot$data$variable, levels=phage_order)
phage_plot$data$CCLabel = factor(phage_plot$data$CCLabel, levels=ccorder)






# 
HGTplot = ggplot(HGT_Info, aes(y=IsolateID, x=variable, color=factor(value), fill=factor(value), starshape=factor(type)))+ geom_star(size=1,size=1) +scale_color_manual(values=colvals)+ scale_starshape_manual(values=starshapeobj)+
  scale_fill_manual(values=fillvals) + theme(axis.text.x= element_text( size=2, angle = 90),  axis.line.y=element_blank(), axis.text.y=element_text(size=2),
                                             legend.position="none", axis.line.x=element_blank(),axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                                             panel.background = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())+scale_x_discrete(position = "top") + coord_fixed(ratio=1)



Phage_plasmid_CC = Phage_plasmid_presence %>% left_join(CCs,by="DORN")
Phage_plasmid_CC$DORN=NULL
Phage_plasmid_CC = Phage_plasmid_CC %>% group_by(CCLabel) %>% summarize_all(max)
saveCCs = (Phage_plasmid_CC$CCLabel)
Phage_plasmid_CC$CCLabel = NULL
TransposedCC_hgts = data.frame(t(Phage_plasmid_CC))
colnames(TransposedCC_hgts) = saveCCs



multiplo_obj = cowplot::plot_grid(phage_plot, plasmid_plot,align = "v",ncol=1)
ggsave(multiplo_obj, file="Documents/DFUData/PhagePlasmidPresence_CC.pdf", height=14.1,width=1.8)


phage_plasmid_pres 
View(phage_plasmid_pres %>% filter(Phage41==1) %>% select(Annotation, Gene))
  
  