library(stringr)
library(stats)
library(tidyr)
library(dplyr)
library(rdist)
library(ggplot2)
library(cluster)

Roary=read.csv("data/Phylogeny2022Data/gene_presence_absence.csv")
Roary[,15:ncol(Roary)] <- if_else(Roary[,15:ncol(Roary)]=="", 0, 1)


View(Roary)
FullData = read.csv("data/staphyloxanthin_paper_data.csv")
BlastHits = read.csv2("data/Phylogeny2022Data/CombinedKinaseTBlastNsFinal.tab", sep="\t",header=F)
colnames(BlastHits) <- c("genomefile", "SubjID", "Something", "PercentID", "LengthHit","SubjStart", "SubjEnd", "QueryCoverage", "Evalue", "Bitscore" )
FullData$DORN = paste0("DORN", FullData$DORN)
# Look for hits at least 80% as long as the reference sak
# And also at least 90% identity and e value <.01
BlastHits$PercentID =sapply(BlastHits$PercentID, function(x) as.numeric(as.character(x)))
BlastHits1 = BlastHits %>% filter(LengthHit > .80*163 & PercentID > 80)
View(BlastHits1)

View(BlastHits)
Resequenced<-c("DORN1197",
               "DORN1340",
               "DORN1471",
               "DORN1602",
               "DORN1646",
               "DORN22",
               "DORN429",
               "DORN801",
               "DORN882",
               "DORN933",
               "DORN929",
               "DORN900",
               "DORN1339",
               "DORN1430",
               "DORN1473",
               "DORN1523",
               "DORN1645",
               "DORN1732",
               "DORN1679",
               "DORN2178")

FullData = FullData %>% mutate(SakUpdated = if_else(DORN %in% BlastHits1$genomefile , "yes", "no"))
FullData = FullData %>% mutate(SakUpdated = if_else( (DORN %in% c("DORN429", "DORN1176")) ,"NoGenome", SakUpdated))
FullData = FullData %>% mutate(ResequencedSince = if_else( (DORN %in% Resequenced) ,"Yes", "No"))

write.csv(FullData, file="data/PaperData_UpdatedSakClassifications.csv")



JustSak = FullData %>% filter(SakUpdated=="yes") %>% filter(!is.na(staphylokinase))

JustSak$Healed = if_else(is.na(JustSak$week_healed), "No", "Yes")
JustSak$HealedBy12 = if_else(is.na(JustSak$week_healed) | JustSak$week_healed>12, "No", "Yes")

t.test(log(staphylokinase)~ Healed, data=JustSak)
t.test(log(staphylokinase)~ HealedBy12, data=JustSak)

FullData = FullData %>% select(-ResequencedSince)
write.csv(FullData, file="data/staphyloxanthin_paper_UpdatedSakClassifications.csv")

RoarySaks = Roary %>% filter(Gene %in% c("sak", "group_2949", "group_2951", "group_2952", "group_2953"))
ContainsPseudogene = Roary %>% filter(Gene =="group_6074") %>% select_if(grepl("DORN", names(.)))
ContainsPseudogene = colnames(ContainsPseudogene %>% select_if(.==1))
mycolsums=colSums(RoarySaks %>% select_if(grepl("DORN", names(.))))


RoarySak = data.frame(mycolsums)
RoarySak$DORN = rownames(RoarySak)
RoarySak %>% filter(DORN %in% ContainsPseudogene)
FullData = FullData %>% left_join(RoarySak, by="DORN")
FullData = FullData %>% mutate(RoaryPredicts = if_else(mycolsums >0, "yes", "no"))
FullData = FullData %>% mutate(RoaryPredicts = if_else(DORN %in% ContainsPseudogene ,"pseudogene", RoaryPredicts))

FullDataGenomesOnly = FullData %>% filter(!is.na(RoaryPredicts))
ggplot(FullDataGenomesOnly, aes(y=log(staphylokinase), x=RoaryPredicts)) + geom_boxplot()  + xlab("Sak gene present")


summary(aov(log(staphylokinase) ~ (RoaryPredicts), data = FullDataGenomesOnly))
