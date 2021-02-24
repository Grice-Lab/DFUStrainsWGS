# Amy Campbell 
# Feb 2021
# Seeing if roary predictions for sak gene presence 
# differ from direct tblastn on 
# or blastn on roary reference 

library(tidyr)
library(ggplot2)
library(dplyr)

# BLAST searches for kinase gene 
blastoutputkinase = data.frame(read.csv("data/CombinedKinaseBlasts.tab", sep='\t', header=F, na.strings=c("", "NA")))
tblastnkinase = data.frame(read.csv("data/CombinedKinaseTBlastNs.tab", sep='\t', header=F))
removed = read.csv("data/ContaminatedIsolates2-12-21.csv", header=T, na.strings=c("", "NA"))
Phenotypedata=read.csv("data/phenotype_variation_11.06.20.csv")

# Roary results -- 125 (14 isolates excluded)
GenePresenceAbsence = read.csv("data/207Isolates/gene_presence_absence.csv", na.strings=c("", "NA"))
GenePresenceAbsence = GenePresenceAbsence %>% filter(Gene=="sak")
kinasepresenceAbsenceDORNs = GenePresenceAbsence[15:ncol(GenePresenceAbsence)]
kinasepresenceAbsenceDORNs = ifelse(is.na(kinasepresenceAbsenceDORNs), 0, 1) 
transposedKinase = data.frame(t(kinasepresenceAbsenceDORNs))
colnames(transposedKinase) = c("sakRoary")
transposedKinase$DORN = rownames(transposedKinase)
Phenotypedata$DORN = paste0("DORN",Phenotypedata$DORN)
PhenotypedataRoary = Phenotypedata %>% left_join(transposedKinase, by="DORN")

# Blast output-- 137
colnames(blastoutputkinase) <- c("genomefile", "SubjID", "Something", "PercentID", "LengthHit","SubjStart", "SubjEnd", "QueryCoverage", "Evalue", "Bitscore" )
blastoutputkinase = blastoutputkinase %>% separate(genomefile, c("Genome", "Cleaned"), extra='drop', sep='_')
blastoutputkinase$Cleaned <- NULL
blastoutputkinase$DORN =blastoutputkinase$Genome
blastoutputkinase = blastoutputkinase %>% filter(! (DORN %in% removed$DORN))
# sak reference was 492 nts; filter out results where length of the hit < 492*.95
blastoutputkinase = blastoutputkinase %>% filter(LengthHit > 492*.95)
blastoutputkinase$DORN

# TblastN for the Uniprot kinase shows the same distribution
colnames(tblastnkinase) <- c("genomefile", "SubjID", "Something", "PercentID", "LengthHit","SubjStart", "SubjEnd", "QueryCoverage", "Evalue", "Bitscore" )
tblastnkinase = tblastnkinase %>% separate(genomefile, c("Genome", "Cleaned"), extra='drop', sep='_')
tblastnkinase$Cleaned <- NULL
tblastnkinase$DORN =tblastnkinase$Genome
tblastnkinase = tblastnkinase %>% filter(! (DORN %in% removed$DORN))
tblastnkinase = tblastnkinase %>% filter(LengthHit > .95*163)

write.table(tblastnkinase$DORN, file="DORNs_Sak_BlastRoary", quote=F, row.names=F, col.names=F)
intersect(tblastnkinase$DORN, transposedKinasepositive$DORN)

PhenotypedataBlast = Phenotypedata %>% left_join(blastoutputkinase, by="DORN")
PhenotypedataBlast = PhenotypedataBlast %>% mutate(KinasePresent = if_else(is.na(PercentID), 0,1))

tblastnkinase$Binary=1
tblastnkinase = tblastnkinase[c("DORN", "Binary")]
length(intersect(tblastnkinase$DORN, blastoutputkinase$DORN))
transposedKinasepositive = transposedKinase %>% filter(sakRoary ==1)

# no 'false positives' in roary's sak assignments
setdiff(transposedKinasepositive$DORN,tblastnkinase$DORN)
# no 'false negatives' either
setdiff(tblastnkinase$DORN,transposedKinasepositive$DORN)

# Testing whether staphylokinase is associated with healing (un-aggregated and aggregated)
PhenotypedataBlastPositives = PhenotypedataBlast %>% filter(KinasePresent==1)
PhenotypedataBlastPositives = PhenotypedataBlastPositives %>% mutate(Healed12 = if_else(WeekHealed > 12, 0, 1))
PhenotypedataBlastPositives = PhenotypedataBlastPositives %>% mutate(HealedYorN = if_else(HealedYorN=="no", 0, 1))

# log transformation doesnt exactly normalize -- should use wilcoxon mann whitney for these 
kinaseplotlog = ggplot(PhenotypedataBlastPositives, aes(x=log(kinase))) + geom_histogram(bins=20, fill="white", color="black")+ ggtitle("Log-transformed Staphylokinase Distribution(sak-positive isolates)")
qqnorm(PhenotypedataBlastPositives$kinase)
qqnorm(log(PhenotypedataBlastPositives$kinase))

# Repeating ARM's association testing but with filter to sak-positive isolates 
PhenotypedataBlastPositivesUnhealed = PhenotypedataBlastPositives %>% filter(HealedYorN==0)
PhenotypedataBlastPositivesHealed = PhenotypedataBlastPositives %>% filter(HealedYorN==1)
PhenotypedataBlastPositivesUnhealed12 = PhenotypedataBlastPositives %>% filter(Healed12==0)
PhenotypedataBlastPositivesHealed12 = PhenotypedataBlastPositives %>% filter(Healed12==1)

t.test(log(kinase) ~ HealedYorN, data=PhenotypedataBlastPositives)
t.test(log(kinase) ~ Healed12, data=PhenotypedataBlastPositives)
wilcox.test(PhenotypedataBlastPositivesUnhealed$kinase, PhenotypedataBlastPositivesHealed$kinase)
wilcox.test(PhenotypedataBlastPositivesUnhealed12$kinase, PhenotypedataBlastPositivesHealed12$kinase)

PhenotypedataBlastPositivesAveraged = PhenotypedataBlastPositives %>% select(c(kinase, Patient, HealedYorN, Healed12)) %>% group_by(Patient) %>% summarize_all(mean)
PhenotypedataBlastPositivesAveragedUnhealed = PhenotypedataBlastPositivesAveraged %>% filter(HealedYorN==0)
PhenotypedataBlastPositivesAveragedHealed = PhenotypedataBlastPositivesAveraged %>% filter(HealedYorN==1)

t.test(log(kinase)~HealedYorN, data=PhenotypedataBlastPositivesAveraged)
wilcox.test(PhenotypedataBlastPositivesAveragedHealed$kinase, PhenotypedataBlastPositivesAveragedUnhealed$kinase)


