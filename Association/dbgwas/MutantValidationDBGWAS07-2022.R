# Amy Campbell
# Lexy, Amelia, and I did xanthin assays on 7 KO mutants for hits from dbgwas 
library(dplyr)
library(ggplot2)

Rep1Xanthin = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/XanthinAssayMutantsJuly2022BioRep1.csv")


Rep1XanthinSummary = Rep1Xanthin %>% group_by(LabelCellPlate) %>% summarise(mean465=round(mean(OD465), 10), mean600=round(mean(OD600), 10 ))
Rep1XanthinSummary$ODratio = round(Rep1XanthinSummary$mean465/Rep1XanthinSummary$mean600, 10)
Rep1XanthinSummary$Norm502A = round(100*((Rep1XanthinSummary$ODratio)/0.22610214), 10)

###################
# Staphyloxanthin
###################

# Replicate 1 
STX1 = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/Validations/STX_1/stx_07.07.22_Values_LexyRep1.csv")
STX1Means = STX1 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX1Means$Ratio = round(STX1Means$mean465/STX1Means$mean600, 10)
STX1502A = STX1Means %>% filter(Strain=="502A") %>% select(Ratio)
STX1Means$Normalized = 100*((STX1Means$Ratio)/STX1502A$Ratio)

STX1Results = STX1Means %>% select(Strain,Normalized)
STX1Results$Replicate=1


# Replicate 2
STX2 = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/Validations/STX_2/stx_07.08.22_Values_LexyRep2.csv")
STX2Means = STX2 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX2Means$Ratio = round(STX2Means$mean465/STX2Means$mean600, 10)
STX2502A = STX2Means %>% filter(Strain=="502A") %>% select(Ratio)
STX2Means$Normalized = 100*((STX2Means$Ratio)/STX2502A$Ratio)
STX2Results = STX2Means %>% select(Strain,Normalized)
STX2Results$Replicate=2

# Replicate 3
STX3 = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/Validations/STX_3/stx_07.12.22_Values_LexyRep3.csv")
STX3Means = STX3 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX3Means$Ratio = round(STX3Means$mean465/STX3Means$mean600, 10)
STX3502A = STX3Means %>% filter(Strain=="502A") %>% select(Ratio)
STX3Means$Normalized = 100*((STX3Means$Ratio)/STX3502A$Ratio)
STX3Results = STX3Means %>% select(Strain, Normalized)
STX3Results$Replicate=3

# Replicate 4
STX4 = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/Validations/STX_4/stx_07.13.22_Values_LexyRep4.csv")
STX4Means = STX4 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX4Means$Ratio = round(STX4Means$mean465/STX4Means$mean600, 10)
STX4502A = STX4Means %>% filter(Strain=="502A") %>% select(Ratio)
STX4Means$Normalized = 100*((STX4Means$Ratio)/STX4502A$Ratio)
STX4Results = STX4Means %>% select(Strain, Normalized)
STX4Results$Replicate = 4
  
# Replicate 5
STX5 = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/Validations/STX_5/stx_07.14.22_Values_LexyRep5.csv")
STX5Means = STX5 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX5Means$Ratio = round(STX5Means$mean465/STX5Means$mean600, 10)
STX5502A = STX5Means %>% filter(Strain=="502A") %>% select(Ratio)
STX5Means$Normalized = 100*((STX5Means$Ratio)/STX5502A$Ratio)
STX5Results = STX5Means %>% select(Strain, Normalized)
STX5Results$Replicate = 5


CombinedReps = do.call("rbind", list(STX1Results, STX2Results, STX3Results, STX4Results, STX5Results))
CombinedReps = CombinedReps %>% filter(! (Strain %in% c("SA113", "SHE", "502A")))
XanthinKOPlot = ggplot(CombinedReps, aes(x=Strain, y=log10(Normalized))) + geom_boxplot(fill="darkgoldenrod") +  ggpubr::stat_compare_means(method = "t.test",ref.group="LAC") + xlab("Strain #") + ylab("Log10(Staphyloxanthin as % 502A)")
orderXaxis = append("LAC", setdiff(unique(CombinedReps$Strain), c("LAC", "SA113", "SHE", "502A")))

orderXaxis  = append(orderXaxis, c("SA113", "SHE", "502A"))
XanthinKOPlot$data$Strain = factor(XanthinKOPlot$data$Strain, levels=orderXaxis)
XanthinKOPlot + theme_classic()


##########
# Biofilm
##########

Biofilm1 = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/Validations/BIOFILM_1/biofilm_07.13.22_Values_LexyRep1.csv")
Biofilm1Means = Biofilm1 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm1Means$Ratio = round(Biofilm1Means$mean570/Biofilm1Means$mean600, 10)
biofilm1_502A = Biofilm1Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm1Means$Normalized = 100*((Biofilm1Means$Ratio)/biofilm1_502A$Ratio)
bfilm1 = Biofilm1Means %>% select(Strain, Normalized)



Biofilm2 = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/Validations/BIOFILM_2/biofilm_07.14.22_Values_LexyRep2.csv")

Biofilm2Means = Biofilm2 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm2Means$Ratio = round(Biofilm2Means$mean570/Biofilm2Means$mean600, 10)
Biofilm2_502A = Biofilm2Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm2Means$Normalized = 100*((Biofilm2Means$Ratio)/Biofilm2_502A$Ratio)
bfilm2 = Biofilm2Means %>% select(Strain, Normalized)

BiofilmCombined = rbind(bfilm1, bfilm2)
BiofilmCombined = BiofilmCombined %>% filter(! (Strain %in% c("SA113", "SHE", "502A")))

ggplot(BiofilmCombined, aes(x=Strain, y=log10(Normalized))) + geom_boxplot()
