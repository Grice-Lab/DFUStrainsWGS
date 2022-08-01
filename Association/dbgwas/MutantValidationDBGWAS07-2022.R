# Amy Campbell
# Lexy, Amelia, and I did xanthin assays on 7 KO mutants for hits from dbgwas 
library(dplyr)
library(ggplot2)

pathstart = "~/Documents/DFUData/Validations/"

###################
# Staphyloxanthin
###################

# Replicate 1 
STX1 = read.csv(paste0(pathstart,"STX_1/stx_07.07.22_Values_LexyRep1.csv"))
STX1Means = STX1 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX1Means$Ratio = round(STX1Means$mean465/STX1Means$mean600, 10)
STX1502A = STX1Means %>% filter(Strain=="502A") %>% select(Ratio)
STX1Means$Normalized = 100*((STX1Means$Ratio)/STX1502A$Ratio)

STX_1Results = STX1Means %>% select(Strain,Normalized)
STX_1Results$Replicate=1


# Replicate 2
STX2 = read.csv(paste0(pathstart,"STX_2/stx_07.08.22_Values_LexyRep2.csv"))
STX2Means = STX2 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX2Means$Ratio = round(STX2Means$mean465/STX2Means$mean600, 10)
STX2502A = STX2Means %>% filter(Strain=="502A") %>% select(Ratio)
STX2Means$Normalized = 100*((STX2Means$Ratio)/STX2502A$Ratio)
STX_2Results = STX2Means %>% select(Strain,Normalized)
STX_2Results$Replicate=2

# Replicate 3
STX3 = read.csv(paste0(pathstart,"STX_3/stx_07.12.22_Values_LexyRep3.csv"))
STX3Means = STX3 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX3Means$Ratio = round(STX3Means$mean465/STX3Means$mean600, 10)
STX3502A = STX3Means %>% filter(Strain=="502A") %>% select(Ratio)
STX3Means$Normalized = 100*((STX3Means$Ratio)/STX3502A$Ratio)
STX_3Results = STX3Means %>% select(Strain, Normalized)
STX_3Results$Replicate=3

# Replicate 4
STX4 = read.csv(paste0(pathstart, "STX_4/stx_07.13.22_Values_LexyRep4.csv"))
STX4Means = STX4 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX4Means$Ratio = round(STX4Means$mean465/STX4Means$mean600, 10)
STX4502A = STX4Means %>% filter(Strain=="502A") %>% select(Ratio)
STX4Means$Normalized = 100*((STX4Means$Ratio)/STX4502A$Ratio)
STX_4Results = STX4Means %>% select(Strain, Normalized)
STX_4Results$Replicate = 4
  
# Replicate 5
STX5 = read.csv(paste0(pathstart, "STX_5/stx_07.14.22_Values_LexyRep5.csv"))
STX5Means = STX5 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX5Means$Ratio = round(STX5Means$mean465/STX5Means$mean600, 10)
STX5502A = STX5Means %>% filter(Strain=="502A") %>% select(Ratio)
STX5Means$Normalized = 100*((STX5Means$Ratio)/STX5502A$Ratio)
STX_5Results = STX5Means %>% select(Strain, Normalized)
STX_5Results$Replicate = 5

# Replicate 6
STX_6 = read.csv(paste0(pathstart, "STX_6/stx_07.19.22_Values_LexyRep6.csv"))
STX_6$Strain = if_else(STX_6$Strain=="She", "SHE", STX_6$Strain)
STX_6Means = STX_6 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX_6Means$Ratio = round(STX_6Means$mean465/STX_6Means$mean600, 10)
STX_6502A = STX_6Means %>% filter(Strain=="502A") %>% select(Ratio)
STX_6Means$Normalized = 100*((STX_6Means$Ratio)/STX_6502A$Ratio)
STX_6Results = STX_6Means %>% select(Strain, Normalized)
STX_6Results$Replicate = 6

# Replicate 7
STX_7 = read.csv(paste0(pathstart, "STX_7/stx_07.20.22_Values_LexyRep7.csv"))
STX_7$Strain = if_else(STX_7$Strain=="She", "SHE", STX_7$Strain)

STX_7Means = STX_7 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX_7Means$Ratio = round(STX_7Means$mean465/STX_7Means$mean600, 10)
STX_7502A = STX_7Means %>% filter(Strain=="502A") %>% select(Ratio)
STX_7Means$Normalized = 100*((STX_7Means$Ratio)/STX_7502A$Ratio)
STX_7Results = STX_7Means %>% select(Strain, Normalized)
STX_7Results$Replicate = 7

# Replicate 8
STX_8 = read.csv(paste0(pathstart, "STX_8/stx_7.21.22_Values_LexyRep8.csv"))

STX_8$Strain = if_else(STX_8$Strain=="She", "SHE", STX_8$Strain)

STX_8Means = STX_8 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX_8Means$Ratio = round(STX_8Means$mean465/STX_8Means$mean600, 10)
STX_8502A = STX_8Means %>% filter(Strain=="502A") %>% select(Ratio)
STX_8Means$Normalized = 100*((STX_8Means$Ratio)/STX_8502A$Ratio)
STX_8Results = STX_8Means %>% select(Strain, Normalized)
STX_8Results$Replicate = 8


# Replicate 9
STX_9 = read.csv(paste0(pathstart, "STX_9/stx_07.22.22_Values_LexyRep9.csv"))
STX_9$Strain = if_else(STX_9$Strain=="She", "SHE", STX_9$Strain)

STX_9Means = STX_9 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX_9Means$Ratio = round(STX_9Means$mean465/STX_9Means$mean600, 10)
STX_9502A = STX_9Means %>% filter(Strain=="502A") %>% select(Ratio)
STX_9Means$Normalized = 100*((STX_9Means$Ratio)/STX_9502A$Ratio)
STX_9Results = STX_9Means %>% select(Strain, Normalized)
STX_9Results$Replicate = 9

# Replicate 10
STX_10 = read.csv(paste0(pathstart, "STX_10/stx_07.27.22_Values_LexyRep10.csv"))
STX_10$Strain = if_else(STX_10$Strain=="She", "SHE", STX_10$Strain)

STX_10Means = STX_10 %>% group_by(Strain) %>% summarise(mean465= round(mean(OD465), 10), mean600= round(mean(OD600), 10))
STX_10Means$Ratio = round(STX_10Means$mean465/STX_10Means$mean600, 10)
STX_10502A = STX_10Means %>% filter(Strain=="502A") %>% select(Ratio)
STX_10Means$Normalized = 100*((STX_10Means$Ratio)/STX_10502A$Ratio)
STX_10Results = STX_10Means %>% select(Strain, Normalized)
STX_10Results$Replicate = 10





CombinedReps = do.call("rbind", list(STX_1Results, STX_2Results, STX_3Results, STX_4Results, STX_5Results, STX_6Results, STX_7Results, STX_8Results, STX_9Results, STX_10Results))
CombinedReps = CombinedReps %>% filter(! (Strain %in% c("SA113", "SHE", "502A")))
CombinedReps$log_10_xanthin = log10(CombinedReps$Normalized)
XanthinKOPlot = ggplot(CombinedReps, aes(x=Strain, y=log_10_xanthin)) + geom_boxplot(fill="darkgoldenrod") +  ggpubr::stat_compare_means(method = "t.test",ref.group="LAC",label= "p.format") + xlab("Strain #") + ylab("Log10(Staphyloxanthin as % 502A)")
orderXaxis = append("LAC", setdiff(unique(CombinedReps$Strain), c("LAC", "SA113", "SHE", "502A")))

orderXaxis  = append(orderXaxis, c("SA113", "SHE", "502A"))
XanthinKOPlot$data$Strain = factor(XanthinKOPlot$data$Strain, levels=orderXaxis)
XanthinKOPlot = XanthinKOPlot +theme_classic() + theme(plot.title=element_text(hjust=.5, size=14, face="bold")) + ggtitle("Staphyloxanthin Tn KOs")


##########
# Biofilm
##########

Biofilm1 = read.csv(paste0(pathstart, "BIOFILM_1/biofilm_07.13.22_Values_LexyRep1.csv"))
Biofilm1Means = Biofilm1 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm1Means$Ratio = round(Biofilm1Means$mean570/Biofilm1Means$mean600, 10)
biofilm1_502A = Biofilm1Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm1Means$Normalized = 100*((Biofilm1Means$Ratio)/biofilm1_502A$Ratio)
bfilm1 = Biofilm1Means %>% select(Strain, Normalized)
bfilm1$Replicate = 1

Biofilm2 = read.csv(paste0(pathstart, "BIOFILM_2/biofilm_07.14.22_Values_LexyRep2.csv"))
Biofilm2Means = Biofilm2 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm2Means$Ratio = round(Biofilm2Means$mean570/Biofilm2Means$mean600, 10)
Biofilm2_502A = Biofilm2Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm2Means$Normalized = 100*((Biofilm2Means$Ratio)/Biofilm2_502A$Ratio)
bfilm2 = Biofilm2Means %>% select(Strain, Normalized)
bfilm2$Replicate = 2

Biofilm3 = read.csv(paste0(pathstart, "BIOFILM_3/biofilm_07.20.22_Values_LexyRep3.csv"))
Biofilm3Means = Biofilm3 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm3Means = Biofilm3 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm3Means$Ratio = round(Biofilm3Means$mean570/Biofilm3Means$mean600, 10)
Biofilm3_502A = Biofilm3Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm3Means$Normalized = 100*((Biofilm3Means$Ratio)/Biofilm3_502A$Ratio)
bfilm3 = Biofilm3Means %>% select(Strain, Normalized)
bfilm3$Replicate = 3

Biofilm4 = read.csv(paste0(pathstart, "BIOFILM_4/biofilm_07.21.22_Values_LexyRep4.csv"))
Biofilm4$Strain = if_else(Biofilm4$Strain=="She", "SHE", Biofilm4$Strain)
Biofilm4Means = Biofilm4 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm4Means = Biofilm4 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm4Means$Ratio = round(Biofilm4Means$mean570/Biofilm4Means$mean600, 10)
Biofilm4_502A = Biofilm4Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm4Means$Normalized = 100*((Biofilm4Means$Ratio)/Biofilm4_502A$Ratio)
bfilm4 = Biofilm4Means %>% select(Strain, Normalized)
bfilm4$Replicate = 4

Biofilm5 = read.csv(paste0(pathstart, "BIOFILM_5/biofilm_07.22.22_Values_LexyRep5.csv"))
Biofilm5$Strain = if_else(Biofilm5$Strain=="She", "SHE", Biofilm5$Strain)
Biofilm5Means = Biofilm5 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm5Means = Biofilm5 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm5Means$Ratio = round(Biofilm5Means$mean570/Biofilm5Means$mean600, 10)
Biofilm5_502A = Biofilm5Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm5Means$Normalized = 100*((Biofilm5Means$Ratio)/Biofilm5_502A$Ratio)
bfilm5 = Biofilm5Means %>% select(Strain, Normalized)
bfilm5$Replicate = 5


Biofilm6 = read.csv(paste0(pathstart, "BIOFILM_6/biofilm_07.28.2022_Values_LexyRep6.csv"))
Biofilm6$Strain = if_else(Biofilm6$Strain=="She", "SHE", Biofilm6$Strain)
Biofilm6Means = Biofilm6 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm6Means = Biofilm6 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm6Means$Ratio = round(Biofilm6Means$mean570/Biofilm6Means$mean600, 10)
Biofilm6_502A = Biofilm6Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm6Means$Normalized = 100*((Biofilm6Means$Ratio)/Biofilm6_502A$Ratio)
bfilm6 = Biofilm6Means %>% select(Strain, Normalized)
bfilm6$Replicate = 6

Biofilm7 = read.csv(paste0(pathstart, "BIOFILM_7/biofilm_07.29.22_Values_LexyRep7.csv"))
Biofilm7$Strain = if_else(Biofilm7$Strain=="She", "SHE", Biofilm7$Strain)
Biofilm7Means = Biofilm7 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm7Means = Biofilm7 %>% group_by(Strain) %>% summarise(mean570= round(mean(OD570), 10), mean600= round(mean(OD600), 10))
Biofilm7Means$Ratio = round(Biofilm7Means$mean570/Biofilm7Means$mean600, 10)
Biofilm7_502A = Biofilm7Means %>% filter(Strain=="502A") %>% select(Ratio)
Biofilm7Means$Normalized = 100*((Biofilm7Means$Ratio)/Biofilm7_502A$Ratio)
bfilm7 = Biofilm7Means %>% select(Strain, Normalized)
bfilm7$Replicate = 7



BiofilmCombined = rbind(bfilm1, bfilm2, bfilm3, bfilm4, bfilm5, bfilm6, bfilm7)
BiofilmCombined = BiofilmCombined %>% filter(! (Strain %in% c("SA113", "SHE", "502A")))

BiofilmKO = ggplot(BiofilmCombined, aes(x=Strain, y=log10(Normalized))) + geom_boxplot( fill="#89E983") + theme_classic() + ggpubr::stat_compare_means(method = "t.test",ref.group="LAC",label= "p.format") + xlab("Strain #") + ylab("Log10(Biofilm as % 502A)")
BiofilmKO = BiofilmKO + theme_classic() + theme(plot.title=element_text(hjust=.5, size=14, face="bold")) + ggtitle("Biofilm Tn KOs")
BiofilmKO$data$Strain = factor(BiofilmKO$data$Strain, levels=orderXaxis)
gridExtra::grid.arrange(BiofilmKO, XanthinKOPlot)

