# Amy Campbell
# Updated June 2022
# Figuring out how to deal with multiple measures per patient when the 'healing' outcome is not independent within patient
# One idea: aggregate them (average the measures of each phenotype within each patient)
# this could mute a real bio signal in cases where a patient has one super virulent strain and one non-virulent strain
# Alternatively, if our question is whether 'higher of any of these phenotypes correlates with nonhealing' we could take the max 
# phenotypic measure for each patient
# As has been pointed out to us now, the  risk of this is that there are more isolates within patient in nonhealing patients, 
# so we'd potentially be more likely to see outliers in those patients' isolates by chance alone than in healing patients
# But, as ARM points out, there are 113 isolates from healed patients total (though the avg # per patient is 10.48
# vs. 3.619 in unhealed vs. healed) --> something we could try is to subsample isolates at random from every patient from whom theres >=4
# Then, take the random subset, choose the highest, go w that. OR, maybe a better way is just to act as though any given unhealed sample
# is taken from the same dist of #s of isolates that the healed are; so, randomly select between 1-7 (max # isolates per healed patient)
# from each unhealed sample


# Daniel shin at the core has also suggested modeling the problem as a multilevel LMM
# predicting isolates' Staphyloxanthin production 
# based on Healing outcome of the patient(as the "exposure" variable) repeated measures of isolates by patient
# with multiple measures per patient


setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/")
library("ggplot2")
library("gridExtra")
library("lme4")
library("dplyr")

FullData= read.csv("data/staphyloxanthin_paper_UpdatedSakClassifications.csv")

set.seed(19104)

###############################################
# Variance in phenotypes: healed vs. nonhealed
###############################################

FullData$Healed = if_else(is.na(FullData$week_healed ), "Unhealed", "Healed")
CountsIsolatesPatients = FullData %>% group_by(patient) %>% dplyr::count()

CountsIsolatesPatientsTimepoints = FullData %>% group_by(patient, visit) %>% dplyr::count()
colnames(CountsIsolatesPatientsTimepoints) = c("patient","visit", "N_isolates_PatientTimepoint")
colnames(CountsIsolatesPatients) = c("patient", "N_isolates_patient")
FullData = FullData %>% left_join(CountsIsolatesPatientsTimepoints, by=c("patient", "visit"))
FullData = FullData %>% left_join(CountsIsolatesPatients, by="patient")
table( (FullData %>% select(patient, Healed) %>% unique() )$Healed)

hist( (FullData %>% filter(Healed=="Healed"))$N_isolates_patient )
hist( (FullData %>% filter(Healed=="Unhealed"))$N_isolates_patient )

median( (FullData %>% filter(Healed=="Healed"))$N_isolates_patient )
median( (FullData %>% filter(Healed=="Unhealed"))$N_isolates_patient )

max(( (FullData %>% filter(Healed=="Unhealed"))$N_isolates_patient ))

# Random subselection test
# Want to randomly subsample nonhealing isolates get *similar* distribution of # of isolates
# Per patient ; basically, go through each patient. if that patient has 4 or more isolates
# and is non-healing, randomly sample between 1-7 isolates(randomly select the # to sample) 
# from that patient for inclusion 
# then, do the 'max xanthin'
# The only thing is, in this case, are you left more likely to get 'extreme' values 
# in the healing group, since there are 2.5X more 
############################################################################################


# spread of phenotypes within healed, non-healed group
#####################################################
FullDataHealed = FullData %>% filter(Healed=="Healed")
FullDataUnhealed = FullData %>% filter(Healed=="Unhealed")

#  staphyloxanthin
sd(log(FullDataHealed$staphyloxanthin))
sd(log(FullDataUnhealed$staphyloxanthin))
IQR(log(FullDataHealed$staphyloxanthin))
IQR(log(FullDataUnhealed$staphyloxanthin))

max(log(FullDataHealed$staphyloxanthin)) - min(log(FullDataHealed$staphyloxanthin)) 
max(log(FullDataUnhealed$staphyloxanthin)) - min(log(FullDataUnhealed$staphyloxanthin)) 

# biofilm
sd(log(FullDataHealed$biofilm))
sd(log(FullDataUnhealed$biofilm))
IQR(log(FullDataHealed$biofilm))
IQR(log(FullDataUnhealed$biofilm))

max(log(FullDataHealed$biofilm)) - min(log(FullDataHealed$biofilm)) 
max(log(FullDataUnhealed$biofilm)) - min(log(FullDataUnhealed$biofilm)) 



# siderophore  (not normally distributed so can't rly interpret SD)
sd((FullDataHealed$siderophore))
sd((FullDataUnhealed$siderophore))
IQR((FullDataHealed$siderophore))
IQR((FullDataUnhealed$siderophore))

max((FullDataHealed$siderophore)) - min((FullDataHealed$siderophore)) 
max((FullDataUnhealed$siderophore)) - min((FullDataUnhealed$siderophore)) 

# staphylokinase (also bimodal or maybe like negative binomial )
hist(FullData$staphylokinase, breaks=20)
sd(FullDataHealed$staphylokinase, na.rm=T)
sd((FullDataUnhealed$staphylokinase), na.rm=T )
IQR(FullDataHealed$staphylokinase, na.rm=T)
IQR(FullDataUnhealed$staphylokinase, na.rm=T)

max((FullDataHealed$staphylokinase), na.rm=T) - min((FullDataHealed$staphylokinase), na.rm=T) 
max((FullDataUnhealed$staphylokinase), na.rm=T) - min((FullDataUnhealed$staphylokinase), na.rm=T) 


#####################################################################################
# Comparison of means in Healed vs. Nonhealed, subsetting isolates to one per patient 
#####################################################################################


#################
# STAPHYLOXANTHIN
#################

PlotHealing = FullData 
PlotHealing$Healed = if_else(is.na(PlotHealing$week_healed ), "Unhealed", "Healed")
PlotHealing$HealedBy12 = if_else( (is.na(PlotHealing$week_healed) | PlotHealing$week_healed > 12) , "Unhealed", "Healed")
PlotHealingAll = PlotHealing
PlotHealing$Healed = factor(PlotHealing$Healed)


# Aggregating within patient (this kind of
# mutes the signals of the higher) s-xanthin strain types
# which may be simultaneously present with low-xanthin strain types
####################################################################
PlotHealingAggregated = PlotHealing %>% group_by(patient) %>% summarize(patient=patient, Healed=Healed, HealedBy12=HealedBy12, meanXanthin= mean(staphyloxanthin)) %>% unique()


AggregatedByPatient = ggplot(PlotHealingAggregated, aes(x=Healed, y=log(meanXanthin))) + geom_boxplot(fill="darkgoldenrod")+  xlab("Patient DFU Healed") + ylab("Staphyloxanthin Production (normalized log-absorbance)") + ggtitle("Staphyloxanthin Production vs. Healing \n(measures aggregated by patient)") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=20, hjust=.5), axis.text.x=element_text( family="Helvetica", size=14),  axis.text.y=element_text( family="Helvetica", size=12),axis.title=element_text(family="Helvetica",size=18))
Unaggregated = ggplot(PlotHealing, aes(x=Healed, y=log(staphyloxanthin))) + geom_boxplot(fill="darkgoldenrod") + xlab("DFU Healing Outcomes") + ylab("Staphyloxanthin Production (normalized log-absorbance)") + ggtitle("Staphyloxanthin Production vs. Healing (221 Isolates)")  + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=20, hjust=.5),  axis.text.x=element_text( family="Helvetica", size=16),  axis.text.y=element_text( family="Helvetica", size=16), axis.title=element_text(family="Helvetica",size=18))


UnaggregatedT = t.test(log(PlotHealing$staphyloxanthin) ~ PlotHealing$Healed)
UnaggregatedP = paste0("p=", round(UnaggregatedT$p.value, 4))

AggregatedT = t.test(log(PlotHealingAggregated$meanXanthin) ~ PlotHealingAggregated$Healed)
AggregatedP = paste0("p=", round(AggregatedT$p.value, 4))

AggregatedByPatient = AggregatedByPatient + annotate(geom="text", x=1.5, y=4.3, label=AggregatedP, size=8)
Unaggregated = Unaggregated + annotate(geom="text", x=1.5, y=4.3, label=UnaggregatedP, size=8)



# Taking maximum phenotype measure per patient
#####################################################
PlotHealingMaxed = PlotHealing %>% group_by(patient) %>% summarize(patient=patient, Healed=Healed, HealedBy12=HealedBy12, maxXanthin= max(staphyloxanthin), maxKinase=max(staphylokinase), maxBiofilm=max(biofilm), maxSiderophore=max(siderophore)) %>% unique()
hist(log(PlotHealingMaxed$maxXanthin))
hist(log(PlotHealingMaxed$maxKinase))
hist(log(PlotHealingMaxed$maxSiderophore))
hist(log(PlotHealingMaxed$maxBiofilm))


# Hypothesis testing with max phenotype per patient
####################################################
T_Test_XanthinMax = t.test(log(PlotHealingMaxed$maxXanthin) ~ PlotHealingMaxed$Healed) 
T_Test_BiofilmMax = t.test(log(PlotHealingMaxed$maxBiofilm) ~ PlotHealingMaxed$Healed) 
Wilcox_SiderophoreMax = wilcox.test(PlotHealingMaxed$maxSiderophore ~ PlotHealingMaxed$Healed)
Wilcox_KinaseMax = wilcox.test(PlotHealingMaxed$maxKinase ~ PlotHealingMaxed$Healed)


# Plots for max phenotype per patient
#####################################

MaxPerPatientXanthin = ggplot(PlotHealingMaxed, aes(x=Healed, y=log(maxXanthin))) + geom_boxplot(fill="darkgoldenrod")+  xlab("DFU Type") + ylab("log(Normalized Staphyloxanthin Production)") + ggtitle("Staphyloxanthin Production vs. Healing \n(maximum staphyloxanthin measure per patient)") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=20, hjust=.5), axis.text.x=element_text( family="Helvetica", size=16),  axis.text.y=element_text( family="Helvetica", size=16),axis.title=element_text(family="Helvetica",size=18))+
  annotate(geom="text", x=1.5, y=4.6, label=paste0("p=",round(T_Test_XanthinMax$p.value,4)), size=8)
ggsave(MaxPerPatientXanthin, file="data/MaxXanthinHealing.pdf", height=7, width=9)

MaxPerPatientKinase= ggplot(PlotHealingMaxed, aes(x=Healed, y=maxKinase)) + geom_boxplot(fill="#008080")+  xlab("DFU Type") + ylab("Normalized Staphylokinase Production") + ggtitle("Maximum Staphylokinase Production vs. Healing") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=20, hjust=.5), axis.text.x=element_text( family="Helvetica", size=16),  axis.text.y=element_text( family="Helvetica", size=16),axis.title=element_text(family="Helvetica",size=18))+
  annotate(geom="text", x=1.5, y=1.4, label=paste0("Wilcoxon p=",round(Wilcox_KinaseMax$p.value,4)), size=8)
ggsave(MaxPerPatientKinase, file="data/MaxKinaseHealing.pdf", height=7, width=9)

MaxPerPatientBiofilm= ggplot(PlotHealingMaxed, aes(x=Healed, y=log(maxBiofilm))) + geom_boxplot(fill="#008080")+  xlab("DFU Type") + ylab("log(Normalized Biofilm Production)") + ggtitle("Maximum Biofilm Production vs. Healing") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=20, hjust=.5), axis.text.x=element_text( family="Helvetica", size=16),  axis.text.y=element_text( family="Helvetica", size=16),axis.title=element_text(family="Helvetica",size=18))+
  annotate(geom="text", x=1.5, y=5, label=paste0("p=",round(T_Test_BiofilmMax$p.value,4)), size=8)
ggsave(MaxPerPatientBiofilm, file="data/MaxBiofilmHealing.pdf", height=7, width=9)

MaxPerPatientKinase= ggplot(PlotHealingMaxed, aes(x=Healed, y=maxKinase)) + geom_boxplot(fill="#008080")+  xlab("DFU Type") + ylab("Normalized Staphylokinase Production") + ggtitle("Maximum Staphylokinase Production vs. Healing") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=20, hjust=.5), axis.text.x=element_text( family="Helvetica", size=16),  axis.text.y=element_text( family="Helvetica", size=16),axis.title=element_text(family="Helvetica",size=18))+
  annotate(geom="text", x=1.5, y=1.4, label=paste0("Wilcoxon p=",round(Wilcox_KinaseMax$p.value,4)), size=8)
ggsave(MaxPerPatientKinase, file="data/MaxKinaseHealing.pdf", height=7, width=9)


# LMM  
#################
lme4::lmer( staphyloxanthin ~ Healed + (1 | patient), data=FullData, REML = F)
lme4::lmer( siderophore ~ Healed + (1 | patient), data=FullData, REML = F)
lme4::lmer( biofilm ~ Healed + (1 | patient), data=FullData, REML = F)
lme4::lmer( staphylokinase ~ Healed + (1 | patient), data=FullData, REML = F)



# Try random subsetting such that there are around the same # of isolates per patient in healed vs. nonhealed
# Max
##############################################################################################################
set.seed(19104)

FullDataSelect = FullData %>% select(patient, visit, DORN, Healed )

# Maximum of 7 isolates per patient in the healing group, so subsample from that total to get closer to the dist 
countsPerPatientFullDataSelect = FullDataSelect %>% filter(Healed=="Healed") %>% group_by(patient) %>% dplyr::count()
max(countsPerPatientFullDataSelect$n)

SubsampledDF = data.frame()
for(pat in unique(FullData$patient)){
  print(pat)
  subdf = FullDataSelect %>% filter(patient==pat)
  # For the Unhealed patients, randomly subsample between 1 and 7 (where the # is randomly chosen)
  if( unique(subdf$Healed)[1]=="Unhealed"){
    numsample = sample(c(1:7),1)
    if(nrow(subdf) <= numsample){
      kept=(1:nrow(subdf))
    }else{
      kept = sample(1:nrow(subdf), numsample)
      
    }
    subdf = subdf[kept, ]
    
  }
  SubsampledDF = rbind(SubsampledDF, subdf)
}

countsPerPatientSubsampled = SubsampledDF %>% group_by(patient) %>% dplyr::count()
SubsampledDF = SubsampledDF %>% left_join(countsPerPatientSubsampled, by="patient")

# Mean, distribution (visually), median should be similar with this subsetting 
mean( (SubsampledDF %>% filter(Healed=="Healed"))$n )
mean( (SubsampledDF %>% filter(Healed=="Unhealed"))$n )

hist( (SubsampledDF %>% filter(Healed=="Healed"))$n, breaks=10 )
hist(( (SubsampledDF %>% filter(Healed=="Unhealed"))$n ), breaks=10)

# still median slightly higher in Unhealed (4 instead of 3) but mean of healed is 3.6 compared to unhealed's 4 
median( (SubsampledDF %>% filter(Healed=="Healed"))$n)
median(( (SubsampledDF %>% filter(Healed=="Unhealed"))$n))

SubsampledData = FullData %>% filter(DORN %in% SubsampledDF$DORN)

SubsampledMaxed = SubsampledData %>% group_by(patient) %>% summarize(patient=patient, Healed=Healed,  maxXanthin= max(staphyloxanthin), maxKinase=max(staphylokinase), maxBiofilm=max(biofilm), maxSiderophore=max(siderophore)) %>% unique()
t.test(log(SubsampledMaxed$maxXanthin) ~ SubsampledMaxed$Healed)


# RSBU KO PLOT
###############
RSBU = read.csv("data/RsbU_stx_06.16.22.csv", header=T)
RSBU = RSBU %>% filter(name!="IlvD")
RSBU$Staphyloxanthin = 100*(RSBU$reading / 0.22566667)
RSBU = RSBU %>% mutate(Strain = if_else( name=="WT", "wildtype", paste0(paste0("\u0394", name) ))) 
ggplot(RSBU, aes(x=Strain, y=Staphyloxanthin)) + geom_boxplot(fill="darkgoldenrod") + theme_classic() + ylab("Normalized Staphyloxanthin Production") + ylim(0,110)
attach(RSBU)
pairwise.t.test(log(Staphyloxanthin), factor(Strain), data= RSBU, p.adj="BH")

