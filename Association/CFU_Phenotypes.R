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
metadataGardner$HealingPoint = case_when( )
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
# Without averaging across isolates from the same subject 
glmxanthin <- glm(HealedEverBin ~ log10(CFU) + Xanthin, data = CombinedFrame, family=binomial)
glmxanthinweek12 <- glm(HealedBy12Bin ~ log10(CFU) + Xanthin, data = CombinedFrame, family=binomial)

# Averages
CombinedFrameSubsetForAverage = CombinedFrame %>% select(c(-DORN, -HealedBy12, -HealedEver))
Averages = CombinedFrameSubsetForAverage %>% group_by(patient_id) %>% summarise_all(mean)
glmAveragesXanthin = glm(HealedEverBin ~ log10(CFU) + Xanthin, data = Averages, family=binomial)
glmAveragesXanthin12 = glm(HealedBy12Bin ~ log10(CFU) + Xanthin, data = Averages, family=binomial)

# Top CFU strain per patient 
TopCFUStrain = CombinedFrame %>% group_by(patient_id) %>% top_n(1, CFU)
glmTopCFUStrain = glm(HealedEverBin ~ log10(CFU) + Xanthin, data = TopCFUStrain, family=binomial)
glmTopCFUStrain12 = glm(HealedBy12Bin ~ log10(CFU) + Xanthin, data = TopCFUStrain, family=binomial)

# Last visit before healing, amputation, or end of study
CombinedFrameMaxVisit = CombinedFrame %>% group_by(patient_id) %>% top_n(1, visit)
# If subject has multiple
CombinedFrameMaxVisit= CombinedFrameMaxVisit %>% group_by(patient_id) %>% top_n(1, CFU)
glmLastTimepointXanthin = glm(HealedEverBin ~ log10(CFU) + Xanthin, data = CombinedFrameMaxVisit)
summary(glmLastTimepointXanthin)
glmLastTimepointXanthin12 = glm(HealedBy12Bin ~ log10(CFU) + Xanthin, data = CombinedFrameMaxVisit)
summary(glmLastTimepointXanthin12)

# 
# # First visit 
# CombinedFrameVisit1 = CombinedFrame %>% filter(visit==1) %>% group_by(patient_id) %>% top_n(1,CFU)
# glmFirstTimepointXanthin = glm(HealedBy12Bin ~ log10(CFU) + Xanthin, data = CombinedFrameVisit1)
# summary(glmFirstTimepointXanthin)

# Testing kinase
################
glmkinase <- glm(HealedEverBin ~ log10(CFU) + kinase, data = CombinedFrame, family=binomial)
glmkinaseweek12 <- glm(HealedBy12Bin ~ log10(CFU) + kinase, data = CombinedFrame, family=binomial)
glmAveragesKinase =  glm(HealedEverBin ~ log10(CFU) + kinase, data = Averages, family = binomial)
glmAveragesKinase12 =  glm(HealedBy12Bin ~ log10(CFU) + kinase, data = Averages, family = binomial)

summary(glmAveragesKinase)
summary(glmAveragesKinase12)

# Kinase associated with healing overall outcome (healed during study vs. not) but not healing by 12 weeks outcome 
summary(glmkinase)
summary(glmkinaseweek12)


