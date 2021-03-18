# Amy Campbell
# Feb 2021

library(dplyr)
library(ggplot2)
setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")

# Big metadata file 
metadataGardner = read.csv("data/gardner_metadata23DEC1.csv")

# Info about isolates 
IsolatesInfo = read.csv("data/DFU_Staph_aureus_isolates.csv")
IsolatesInfo$DORN = paste0("DORN",IsolatesInfo$Doern.lab.bank.)
IsolatesInfo = IsolatesInfo %>% select(DORN, CFU, patient_id, visit)

# 11/2020 phenotypes
phenotypes11 = read.csv("data/phenotype_variation_11.06.20.csv")
phenotypes11$DORN = paste0("DORN", phenotypes11$DORN)

phenotypes11 = phenotypes11 %>% select(c(DORN,Patient, siderophore, hemolysis, biofilm, kinase, WeekHealed))

# Test metadata compatibility for healing outcome + weeks to heal 
metadataGardner$Patient = metadataGardner$studyid

testMetaDataCompatibility = metadataGardner %>% left_join(phenotypes11[c("Patient", "WeekHealed")])
LK = read.csv("data/metamapLK.csv")
LK = LK %>% select(c(patient_id, healing_time, outcome, visit_healed)) 
colnames(LK) = c("Patient", "healing_time", "outcome", "visit_healed")

metadatagardner_subset = metadataGardner %>% select(c(Patient, wks_to_heal,visit_healed))
metadatagardner_subset$weeks_to_heal = round(metadatagardner_subset$wks_to_heal)

testLK_vs_Gardner = metadatagardner_subset %>% left_join(LK, by="Patient")
testLK_vs_Gardner_McCready = testLK_vs_Gardner %>% left_join(phenotypes11[c("Patient", "WeekHealed")])

# Amelia's method of setting 'weeks to heal' to 50 if unhealed 
metadatagardner_subset$weeks_to_heal[is.na(metadatagardner_subset$weeks_to_heal)] <- 50


# Staphyloxanthin assay finished for those two missing strains 
UptoDateXanthin = read.csv("data/staphyloxanthin_averages_1.13.21.csv")
UptoDateXanthin$DORN = paste0("DORN", UptoDateXanthin$X)
UptoDateXanthin = UptoDateXanthin %>% select(c(DORN, Average))
colnames(UptoDateXanthin) = c("DORN", "Xanthin")

# Merge the three together, starting with IsolatesInfo 
CombinedFrame = IsolatesInfo %>% left_join(phenotypes11, by="DORN")
CombinedFrame = CombinedFrame %>% left_join(UptoDateXanthin, by="DORN")
metadatagardner_subset$patient_id = metadatagardner_subset$Patient
CombinedFrame = CombinedFrame %>% left_join(metadatagardner_subset[c("patient_id", "weeks_to_heal")], by="patient_id")
  
# Filter out the strains that are neither in Amelia's Jan 2020 phenotypes data nor in her updated Xanthin data
CombinedFrame = CombinedFrame %>% filter(!(is.na(Patient) & is.na(Xanthin)))
# remove DORN685 (actually s. xylosus) and DORN946(actually s. simulans)
CombinedFrame = CombinedFrame %>% filter( !(DORN %in% c("DORN685", "DORN946")))
CombinedFrame$HealedBy12 = if_else(CombinedFrame$weeks_to_heal <= 12, "Yes", "No")
CombinedFrame$HealedEver = if_else(CombinedFrame$weeks_to_heal == 50, "No", "Yes")
CombinedFrame$HealedEverBin = if_else(CombinedFrame$HealedEver=="No", 0, 1)
CombinedFrame$HealedBy12Bin = if_else(CombinedFrame$HealedBy12=="No", 0, 1)
CombinedFrame$Xanthin = sapply(CombinedFrame$Xanthin , function(x) as.numeric(as.character(x)))

#########################
# TESTING STAPHYLOXANTHIN 
#########################

# Normality 
xanthinplot_raw = ggplot(CombinedFrame, aes(x=Xanthin)) + geom_histogram(bins=40,fill="white", color="black") + ggtitle("Untransformed Staphyloxanthin Distribution")
xanthinplot_log = ggplot(CombinedFrame, aes(x=log(Xanthin))) + geom_histogram(bins=40, fill="white", color="black")+ ggtitle("Log-transformed Staphyloxanthin Distribution")
gridExtra::grid.arrange(xanthinplot_raw, xanthinplot_log)

xanthinNorm = qqnorm(CombinedFrame$Xanthin, main="Un-transformed Xanthin") 
xanthinNormLog = qqnorm(log(CombinedFrame$Xanthin), main="Log-transformed Xanthin")

t.test(log(Xanthin)~HealedEverBin, data=CombinedFrame)
t.test(log(Xanthin)~HealedBy12Bin, data=CombinedFrame)

# Averages
CombinedFrameSubsetForAverage = CombinedFrame %>% select(c(-DORN, -HealedBy12, -HealedEver))
Averages = CombinedFrameSubsetForAverage %>% group_by(patient_id) %>% summarise_all(mean)
qqnorm(Averages$Xanthin, main="Subject-Averaged Staphyloxanthin(untransformed)")

t.test(log(Xanthin)~HealedEverBin, data=Averages)
t.test(log(Xanthin)~HealedBy12Bin, data=Averages)

# Without averaging across isolates from the same subject 
glmxanthin <- glm(HealedEverBin ~ log(CFU) + log(Xanthin), data = CombinedFrame, family=binomial)
summary(glmxanthin)
glmxanthinweek12 <- glm(HealedBy12Bin ~ log(CFU) + log(Xanthin), data = CombinedFrame, family=binomial)
summary(glmxanthinweek12)

glmAveragesXanthin = glm(HealedEverBin ~ log(CFU) + log(Xanthin), data = Averages, family=binomial)
summary(glmAveragesXanthin)
glmAveragesXanthin12 = glm(HealedBy12Bin ~ log(CFU) + log(Xanthin), data = Averages, family=binomial)
summary(glmAveragesXanthin12)


# Top CFU strain per patient 
TopCFUStrain = CombinedFrame %>% group_by(patient_id) %>% top_n(1, CFU)
glmTopCFUStrain = glm(HealedEverBin ~ log(CFU) + Xanthin, data = TopCFUStrain, family=binomial)
glmTopCFUStrain12 = glm(HealedBy12Bin ~ log(CFU) + Xanthin, data = TopCFUStrain, family=binomial)




# Last visit before healing, amputation, or end of study
CombinedFrameMaxVisit = CombinedFrame %>% group_by(patient_id) %>% top_n(1, visit)
CombinedFrameMaxVisit= CombinedFrameMaxVisit %>% group_by(patient_id) %>% top_n(1, CFU)
glmLastTimepointXanthin = glm(HealedEverBin ~ log(CFU) + log(Xanthin), data = CombinedFrameMaxVisit)
summary(glmLastTimepointXanthin)
glmLastTimepointXanthin12 = glm(HealedBy12Bin ~ log(CFU) + log(Xanthin), data = CombinedFrameMaxVisit)
summary(glmLastTimepointXanthin12)

# 
# # First visit 
# CombinedFrameVisit1 = CombinedFrame %>% filter(visit==1) %>% group_by(patient_id) %>% top_n(1,CFU)
# glmFirstTimepointXanthin = glm(HealedBy12Bin ~ log(CFU) + Xanthin, data = CombinedFrameVisit1)
# summary(glmFirstTimepointXanthin)

# Testing kinase
################

# Normality
# woooooof idk if kinase can be treated as a normally distributed variable ever 
qqnorm(CombinedFrame$kinase)
qqnorm(log(CombinedFrame$kinase))


glmkinase <- glm(HealedEverBin ~ log(CFU) + kinase, data = CombinedFrame, family=binomial)
glmkinaseweek12 <- glm(HealedBy12Bin ~ log(CFU) + kinase, data = CombinedFrame, family=binomial)
glmAveragesKinase =  glm(HealedEverBin ~ log(CFU) + kinase, data = Averages, family = binomial)
glmAveragesKinase12 =  glm(HealedBy12Bin ~ log(CFU) + kinase, data = Averages, family = binomial)

summary(glmAveragesKinase)
summary(glmAveragesKinase12)

# Kinase associated with healing overall outcome (healed during study vs. not) but not healing by 12 weeks outcome 
summary(glmkinase)
summary(glmkinaseweek12)

ggplot(data=CombinedFrame, aes(x=log(Xanthin),y=log(CFU), color=HealedEver)) + geom_point() + scale_color_manual(values=c("#56B4E9", "#D55E00"))
ggplot(data=Averages, aes(x=log(Xanthin),y=log(CFU), color=as.factor(HealedEverBin))) + geom_point() + scale_color_manual(values=c("#56B4E9", "#D55E00"))

CombinedFrameDORNsOrdered = CombinedFrame %>% group_by(patient_id) %>% summarize(meanXanthin=mean(Xanthin)) 

boxplotXanthin=ggplot(data=CombinedFrame, aes(x=as.factor(HealedEverBin), y=log(Xanthin))) + geom_boxplot() 
boxplotXanthin
XanthinBox = ggplot(Averages, aes(x=as.factor(HealedEverBin), y=log(Xanthin))) + geom_boxplot()
