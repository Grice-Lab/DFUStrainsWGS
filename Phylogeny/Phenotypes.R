# Amy Campbell
# 09/29/2020 
# Staphyloxanthin phenotype and healing
# Amelia's found an association between 'healing outcome for the subject a strain came from' 
# and 'that strain's staphyloxanthin phenotype'
# Since there's inherently a bias in treating each strain independently (since subjects with ), 
# Here, I check the association when you average the staphyoxanthin phenotype across all strains present in each subject
# prior to testing for an association (log-transformed to better approximate a normal dist. )

library(dplyr)
library(ggplot2)

# 685 and 946 have been discovered to not be S. aureus (xylosis & simulans)
Biofilm_Xanthin = read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/Phenotypes_01.15.20.csv")
Biofilm_Xanthin = Biofilm_Xanthin %>% filter(DORN!=685 & DORN!=946)

Average_BioFilmXanthin = Biofilm_Xanthin[c("Patient", "WeekHealed","BiofilmPhenotype","XanthinPhenotype")] %>% group_by(Patient) %>% summarise_all(mean)
Average_BioFilmXanthin$HealedBy12 = if_else(Average_BioFilmXanthin$WeekHealed > 12, FALSE, TRUE)

# Transformations necessary
qqnorm(Average_BioFilmXanthin$XanthinPhenotype)
qqnorm(log10(Average_BioFilmXanthin$XanthinPhenotype))
Average_BioFilmXanthin$Log10Xanthin = log10(Average_BioFilmXanthin$XanthinPhenotype)

# Lindsay's paper used the healed by 12 weeks or didn't heal until >12 weeks criteria 
# test Xanthin association with healing outcome under these criteria 
t.test(Log10Xanthin~HealedBy12, data=Average_BioFilmXanthin)
ggplot(Average_BioFilmXanthin, aes(HealedBy12, Log10Xanthin)) + geom_boxplot()

# Amelia has thus far been using the 'healed ever' criteria
# When you average the strains' staphyloxanthin 
Average_BioFilmXanthin$HealedEver= if_else(Average_BioFilmXanthin$WeekHealed >=50, FALSE, TRUE)
t.test(Log10Xanthin~HealedEver, data=Average_BioFilmXanthin)
ggplot(Average_BioFilmXanthin, aes(HealedEver, Log10Xanthin)) + geom_boxplot()

# Figuring out which DORNs are left out from amelia's phenotyping analyses and which are left out from my 
BigFrameDORNS = sapply(bigframe$DORN, function(x) as.numeric(as.character(stringr::str_remove(x,"DORN"))))
setdiff(BigFrameDORNS, Biofilm_Xanthin$DORN)
setdiff( Biofilm_Xanthin$DORN, BigFrameDORNS)

# Subjects with multiple S. aureus lineages (>250 snps in core genome)
#Looking at the staphyloxanthin phenotypes for some subjects who've been shown to have multiple divergent lineages at once
#Specifically, I wonder if 

numgroups = read.csv("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/PhylogeneticGroups_By_Patient.csv")
numgroups_nonzero = numgroups %>% filter(n>1)
PatientOutcomesMultigroups = Average_BioFilmXanthin %>%filter (Patient %in% (numgroups_nonzero$patient_id))
Biofilm_Xanthin %>% filter(Patient %in% c("141", "176", "111", "120", "124", "159", "191"))


groupAssigns = read.csv("/Users/amycampbell/Desktop/GriceLabGit/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/PhylogeneticGroups_By_Patient.csv")
