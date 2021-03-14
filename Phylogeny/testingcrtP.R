# Amy Campbell
# March 2021
# Trying to figure out what's up with this missing chunk (length 174) at the end of some of the crtP (in the crtP gene alignment from PRANK)
# Protein sequence is the same up to ...EDIEK, then has this weird IIVLIVVQYMVL sequence and then ends (as opposed to going for the 70 AA the reference does)

library(seqinr)
library(dplyr)
library(ggplot2)
genomeDORN2127=read.fasta("/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning_221GenomesArchive/RoaryOutput/RoaryResults/ProkkaResults/DORN2127/DORN2127.fna")

# # crtP is located:
# 79145:77826 in feature 1
feature1 = (genomeDORN2127$'1')
crtPsequence= feature1[77826:79145]

# There's a 174-length chunk 'missing' from the crtP 
crtPsequencePlus174 = feature1[(77826-174):79145]


# Genome DORN568 has crtP on feature 10 from 41110:39617
# 41110	39617 DORN568 CrtP Feature 10
genomeDORN568 = read.fasta("/Volumes/QuarantineBackup/DFU_Saureus_PostCleaning_221GenomesArchive/RoaryOutput/RoaryResults/ProkkaResults/DORN568/DORN568.fna")
feature10 = genomeDORN568$'10'
crtPsequenceFull = feature10[39617:41110]
length(crtPsequenceFull)
length(crtPsequencePlus164)

crtPsequenceFull[1:173]
crtPsequencePlus164[2:174]

df_compare = data.frame(Dorn568= crtPsequenceFull,Dorn2127 = crtPsequencePlus164)
df_compare$equal = if_else(df_compare$Dorn568==df_compare$Dorn2127,"True", "False")

View(df_compare)
# bases mostly equivalent (equal=True) from rows 1494 to row 213, then at row 212 everything goes awry
# In row 212 of this dataframe, there's a 't' missing in DORN2127 which causes a frameshift such that instead of E(219:217), K(216:214), N(213:211), Y(210:208)
# it's E(219:217) K(216:214) I I V L I V V Q Y M V L STOP


# all of these are in the CC5 lineage-- are they all from the same subject? 
crtPMutants = c("DORN2127", "DORN2004",
"DORN1869", "DORN1902",
"DORN2034", "DORN2205",
"DORN2166", "DORN2139",
"DORN1952", "DORN807",
"DORN1849", "DORN2075",
"DORN1176", "DORN1829",
"DORN1968", "DORN1844")

Isolate_Info = read.csv("/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/DataInput_DFUStrainWGS_scripts/data/DFU_Staph_aureus_isolates.csv")

Isolate_Info$DORN = paste0("DORN", Isolate_Info$Doern.lab.bank.)
Isolate_Info_Mutants = Isolate_Info %>% filter(DORN %in% crtPMutants)
View(Isolate_Info_Mutants)
patients = unique(Isolate_Info_Mutants$patient_id)

gardnermetadata = read.csv("/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/DataInput_DFUStrainWGS_scripts/data/gardner_metadata23DEC1.csv")
gardnermetadata = gardnermetadata %>% filter(studyid %in% patients)

# 2 of them (136 and 145) healed
# 2 of them didn't (186 still wasn't healed by the end; 197 was amputated)

# CC5 clade
CC5s = read.csv("/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/DataInput_DFUStrainWGS_scripts/data/207Isolates/NodeChildren_Collapsed2000_1.csv")
CC5s = CC5s %>% mutate(MutantCrtP = if_else(DORN %in% crtPMutants, "True", "False"))

# Staphyloxanthin measures 
xanthins = read.csv("/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/DataInput_DFUStrainWGS_scripts/data/staphyloxanthin_averages_1.13.21.csv")
xanthins$DORN = paste0("DORN", xanthins$X)
xanthins = xanthins[c("DORN", "Average")]
CC5s = CC5s %>% left_join(xanthins, by="DORN")

# 
CC5s$Average = as.numeric(as.character(CC5s$Average))
CC5s$MutantCrtP = factor(CC5s$MutantCrtP)

# NA rows are internal nodes & references
test_xanthin_CC5_crtp = t.test(log(as.numeric(as.character(Average))) ~ factor(MutantCrtP), data=CC5s)

ggplot(CC5s, aes(x=MutantCrtP, y=log(Average))) + geom_boxplot( fill="darkgoldenrod1") + xlab("crtP Reading Frame Mutation") + ylab("Log Stapyloxanthin OD") + ggtitle("Staphyloxanthin Production: crtP mutation (CC5 Lineage Only)") + annotate("text", label=paste0("p=", toString(test_xanthin_CC5_crtp$p.value)), y= 4.5, x=1.5, size=5)
