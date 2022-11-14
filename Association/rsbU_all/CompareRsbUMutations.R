library(dplyr)
library(ggplot2)
# Quick and dirty comparison of rsbU NS mutation vs. identical AA sequence to FPRP3757


BlastPBrowserResults = read.csv("/Users/amycampbell/Documents/GriceLabGit/DFUStrainsWGS/Association/rsbU_all/BlastPRsbuPresAbsenceCompare.csv")
Phenotypes= read.csv("/Users/amycampbell/Documents/DFUData/StaphyloxanthinPaperData/staphyloxanthin_paper_phenotypes.csv")

Phenotypes$ID = Phenotypes$DORN


STXvsRsbU = Phenotypes %>% select(ID, staphyloxanthin) %>% left_join(BlastPBrowserResults %>% select(ID,AA_Identity ), by="ID")
STXvsRsbU = STXvsRsbU %>% filter(!is.na(AA_Identity))
ggplot(STXvsRsbU, aes(x=AA_Identity, y=log(staphyloxanthin))) + geom_boxplot() + ggpubr::stat_compare_means(method="t.test")


MismatchPositions = read.csv("/Users/amycampbell/Documents/GriceLabGit/DFUStrainsWGS/Association/rsbU_all/RsbUMismatches.csv")

CCs = read.csv("Documents/DFUData/CCMapPlotting.csv")

CCs$Genome = CCs$DORN

MismatchPositions = MismatchPositions %>% left_join(CCs, by="Genome")


ggplot(MismatchPositions, aes(x=NonmatchPos, fill=CCLabel)) + geom_bar() + theme_classic() + scale_fill_brewer(palette="Dark2")


Mismatch_in_Phosphatase= MismatchPositions  %>% filter(NonmatchPos>= 120 & NonmatchPos<=331)


STXvsRsbU = STXvsRsbU %>% mutate(MismatchPhos = if_else(ID %in% Mismatch_in_Phosphatase$DORN, "yes", "no"))
ggplot(STXvsRsbU, aes(x=MismatchPhos, y=log(staphyloxanthin))) + geom_boxplot(fill="goldenrod") + ggpubr::stat_compare_means(method="t.test")  + theme_classic() + xlab("Mismatched AA to LAC in Phosphatase Domain (120-331)") + ylab("log(staphyloxanthin) as % 502A")


AllIdenticals = BlastPBrowserResults %>% filter(AA_Identity=="Identical")

AllIdenticals$Genome= AllIdenticals$ID

AllIdenticals = AllIdenticals %>% select(Genome) %>% left_join(CCs,by="Genome")

View(AllIdenticals)


MismatchesCC5 = MismatchPositions %>% filter(CCLabel=="CC5")
table(MismatchesCC5$Genome)

View(MismatchesCC5)
