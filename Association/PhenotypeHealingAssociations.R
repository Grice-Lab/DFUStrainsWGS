# Amy Campbell
# 2022
# Figuring out how to deal with multiple measures per patient when the 'healing' outcome is not independent within patient
# One idea: aggregate them (average the measures of each phenotype within each patient)
# this could mute a real bio signal in cases where a patient has one super virulent strain and one non-virulent strain
# Alternatively, if our quesiton is whether 'higher of any of these phenotypes correlates with nonhealing' we could take the max 
# phenotypic measure for each patient

#################
# STAPHYLOXANTHIN
#################

PlotHealing = FullData 
PlotHealing$Healed = if_else(is.na(PlotHealing$week_healed ), "Unhealed", "Healed")
PlotHealing$HealedBy12 = if_else( (is.na(PlotHealing$week_healed) | PlotHealing$week_healed > 12) , "Unhealed", "Healed")
PlotHealingAll = PlotHealing
PlotHealing$Healed = factor(PlotHealing$Healed)
ggplot(PlotHealing, aes(x=Healed, y=log(staphyloxanthin))) + geom_boxplot()

PlotHealing = PlotHealing %>% select(staphyloxanthin, Healed,HealedBy12, patient)

# Aggregating within patient
############################
PlotHealingAggregated = PlotHealing %>% group_by(patient) %>% summarize(patient=patient, Healed=Healed, HealedBy12=HealedBy12, meanXanthin= mean(staphyloxanthin)) %>% unique()
hist(log(PlotHealingAggregated$meanXanthin))


AggregatedByPatient = ggplot(PlotHealingAggregated, aes(x=Healed, y=log(meanXanthin))) + geom_boxplot(fill="darkgoldenrod")+  xlab("Patient DFU Healed") + ylab("Staphyloxanthin Production (mormalized log-absorbance)") + ggtitle("Staphyloxanthin Production vs. Healing \n(measures aggregated by patient)") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=24, hjust=.5), axis.text.x=element_text( family="Helvetica", size=12),  axis.text.y=element_text( family="Helvetica", size=12),axis.title=element_text(family="Helvetica",size=18))
Unaggregated = ggplot(PlotHealing, aes(x=Healed, y=log(staphyloxanthin))) + geom_boxplot(fill="darkgoldenrod") + xlab("DFU Healing Outcome at 26 Weeks") + ylab("Staphyloxanthin Production (mormalized log-absorbance)") + ggtitle("Staphyloxanthin Production vs. Healing")  + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=24, hjust=.5),  axis.text.x=element_text( family="Helvetica", size=14),  axis.text.y=element_text( family="Helvetica", size=14), axis.title=element_text(family="Helvetica",size=18))

UnaggregatedT = t.test(log(PlotHealing$staphyloxanthin) ~ PlotHealing$Healed)
UnaggregatedP = paste0("p=", round(UnaggregatedT$p.value, 4))

AggregatedT = t.test(log(PlotHealingAggregated$meanXanthin) ~ PlotHealingAggregated$Healed)
AggregatedP = paste0("p=", round(AggregatedT$p.value, 4))

AggregatedByPatient = AggregatedByPatient + annotate(geom="text", x=1.5, y=4.3, label=AggregatedP, size=8)
Unaggregated = Unaggregated + annotate(geom="text", x=1.5, y=4.3, label=UnaggregatedP, size=8)

gridExtra::grid.arrange(Unaggregated, AggregatedByPatient, ncol=2)


# Taking maximum staphyloxanthin measure per patient
#####################################################
PlotHealingMaxed = PlotHealing %>% group_by(patient) %>% summarize(patient=patient, Healed=Healed, HealedBy12=HealedBy12, maxXanthin= max(staphyloxanthin), maxKinase=max(staphylokinase)) %>% unique()
hist(log(PlotHealingMaxed$maxXanthin))

MaxPerPatient = ggplot(PlotHealingMaxed, aes(x=Healed, y=log(maxXanthin))) + geom_boxplot(fill="darkgoldenrod")+  xlab("Patient DFU Healed") + ylab("Staphyloxanthin Production (mormalized log-absorbance)") + ggtitle("Staphyloxanthin Production vs. Healing \n(maximum staphyloxanthin measure per patient)") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=24, hjust=.5), axis.text.x=element_text( family="Helvetica", size=12),  axis.text.y=element_text( family="Helvetica", size=12),axis.title=element_text(family="Helvetica",size=18))+
  annotate(geom="text", x=1.5, y=4.3, label=paste0("p=",round(MaxedT$p.value,3)), size=8)

MaxedT = t.test(log(PlotHealingMaxed$maxXanthin) ~ PlotHealingMaxed$Healed) 

MaxedT12 = t.test(log(PlotHealingMaxed$maxXanthin) ~ PlotHealingMaxed$HealedBy12) 
MaxPerPatient12 = ggplot(PlotHealingMaxed, aes(x=HealedBy12, y=log(maxXanthin))) + geom_boxplot(fill="darkgoldenrod")+  xlab("Patient DFU Healed") + ylab("Staphyloxanthin Production (mormalized log-absorbance)") + ggtitle("Staphyloxanthin Production vs. Healing by 12 Weeks \n(maximum staphyloxanthin measure per patient)") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=24, hjust=.5), axis.text.x=element_text( family="Helvetica", size=12),  axis.text.y=element_text( family="Helvetica", size=12),axis.title=element_text(family="Helvetica",size=18))+
  annotate(geom="text", x=1.5, y=4.3, label=paste0("p=",round(MaxedT12$p.value,3)), size=8)

################
# STAPHYLOKINASE
################

UpdatedSak = read.csv("data/staphyloxanthin_paper_UpdatedSakClassifications.csv")
sakYes = UpdatedSak %>% filter(SakUpdated=="yes")
sakYes$Healed = if_else(is.na(sakYes$week_healed ), "Unhealed", "Healed")
sakYes$HealedBy12 = if_else(is.na(sakYes$week_healed ) | sakYes$week_healed >12, "Unhealed", "Healed")

# Unaggregated
t.test(log(sakYes$staphylokinase) ~ sakYes$Healed)
t.test(log(sakYes$staphylokinase) ~ sakYes$HealedBy12)

sakYesMax = sakYes %>% select(staphylokinase, patient, Healed, HealedBy12) %>% group_by(patient) %>% summarise(Healed=Healed, HealedBy12 = HealedBy12, maxKinase=max(staphylokinase))
MaxedSakT = t.test(log(sakYesMax$maxKinase) ~ sakYesMax$Healed)
t.test(log(sakYesMax$maxKinase) ~ sakYesMax$HealedBy12)
ggplot(sakYesMax, aes(x=Healed, y=log(maxKinase))) + geom_boxplot(fill="#008080")+  xlab("Patient DFU Healed") + ylab("log(Staphylokinase Production)") + ggtitle("Maximum Staphylokinase Production vs. Healing \n(sak-positive genomes)") + theme_classic() +
  theme(plot.title=element_text(face="bold", family="Helvetica", size=24, hjust=.5), axis.text.x=element_text( family="Helvetica", size=12),  axis.text.y=element_text( family="Helvetica", size=12),axis.title=element_text(family="Helvetica",size=18))+
  annotate(geom="text", x=1.5, y=.8, label=paste0("p=",round(MaxedSakT$p.value,3)), size=8)



FullData$Healed = if_else(is.na(FullData$week_healed), "Unhealed", "Healed")
FullDataMaxPhenotypes  = FullData %>% select(Healed, patient, staphylokinase, staphyloxanthin) %>% group_by(patient) %>% summarise(MaxStaphylokinase = max(staphylokinase), MaxStaphyloxanthin=max(staphyloxanthin), Healed=Healed, patient=patient)

hist((log((FullDataMaxPhenotypes$MaxStaphylokinase))))
FullDataMaxPhenotypes$MaxStaphylokinaseLog = log(FullDataMaxPhenotypes$MaxStaphylokinase)
FullDataMaxPhenotypes$MaxStaphyloxanthinLog = log(FullDataMaxPhenotypes$MaxStaphyloxanthin ) 
FullDataMaxPhenotypes$MaxStaphylokinaseLog = sapply(FullDataMaxPhenotypes$MaxStaphylokinaseLog, function(x) as.numeric(as.character(x)))
FullDataMaxPhenotypes$MaxStaphyloxanthinLog = sapply(FullDataMaxPhenotypes$MaxStaphyloxanthinLog, function(x) as.numeric(as.character(x)))
FullDataMaxPhenotypes$Healed = factor(FullDataMaxPhenotypes$Healed)


summary(glm(Healed ~   MaxStaphylokinaseLog ,data=FullDataMaxPhenotypes, family="binomial"))

##########
# BIOFILM
##########

PlotHealingBiofilm = FullData 
PlotHealingBiofilm$Healed = if_else(is.na(PlotHealingBiofilm$week_healed ), "Unhealed", "Healed")

PlotHealingBiofilmMaxed = PlotHealingBiofilm %>% select(biofilm, Healed, patient) %>% group_by(patient) %>% summarise(maxBiofilm = max(biofilm), Healed=Healed, patient=patient)
hist(log(PlotHealingBiofilmMaxed$maxBiofilm))
t.test(log(PlotHealingBiofilmMaxed$maxBiofilm) ~ PlotHealingBiofilmMaxed$Healed)

###############
# SIDEROPHORE
################
PlotHealingSid = FullData 
PlotHealingSid$Healed = if_else(is.na(PlotHealingSid$week_healed ), "Unhealed", "Healed")

PlotHealingSiderophoreMaxed = PlotHealingSid %>% select(siderophore, Healed, patient) %>% group_by(patient) %>% summarise(maxSiderophore = max(siderophore), Healed=Healed, patient=patient)
t.test((PlotHealingSiderophoreMaxed$maxSiderophore) ~ PlotHealingSiderophoreMaxed$Healed)
