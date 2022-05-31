library(stringr)
library(stats)
library(tidyr)
library(dplyr)
library(rdist)
library(ggplot2)
library(cluster)

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

