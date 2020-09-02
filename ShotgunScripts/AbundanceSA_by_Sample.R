# Amy Campbell
# 04/2020
# Merging amelia's phenotype dataframes 
library("dplyr")
library("stringr")

# Read in the table that has abundance of each organism in each 
MetaPhlanAbundancetable=read.table('/Users/amycampbell/Desktop/pooled_fulltable.tsv')
MetaPhlanAbundancetable_Staphylococcus=data.frame(t(MetaPhlanAbundancetable))
MetaPhlanAbundancetable_Staphylococcus = MetaPhlanAbundancetable_Staphylococcus["Staphylococcus_aureus"]
MetaPhlanAbundancetable_Staphylococcus$Sample=rownames(MetaPhlanAbundancetable_Staphylococcus)

# Some messy as hell functions applied to the sample names in this abundance table to extract subject and timepoint
MetaPhlanAbundancetable_Staphylococcus$Sample = sapply(MetaPhlanAbundancetable_Staphylococcus$Sample, function(x) stringr::str_replace(x, "filtered_sorted_", ""))
MetaPhlanAbundancetable_Staphylococcus$Timepoint = sapply(MetaPhlanAbundancetable_Staphylococcus$Sample,function(x) as.numeric((strsplit(x, '\\.'))[[1]][2]))
MetaPhlanAbundancetable_Staphylococcus$Subject = sapply(MetaPhlanAbundancetable_Staphylococcus$Sample,function(x) as.numeric((strsplit(x, '\\.'))[[1]][1]))

# Check that all 195 metagenomic samples are represented 
dim(MetaPhlanAbundancetable_Staphylococcus)

MetaPhlanAbundancetable_Staphylococcus$subject_timepoint = paste(MetaPhlanAbundancetable_Staphylococcus$Subject, MetaPhlanAbundancetable_Staphylococcus$Timepoint, sep="_")
write.csv(MetaPhlanAbundancetable_Staphylococcus, "Staphylococcus_by_Metagenome.csv")

# The DOERN biobank isolates (232 of them-- before the cleanup)
StaphIsolateDORNs=read.csv("DFU_Staph_aureus_isolates.csv")
StaphIsolateDORNs$visit_1 = StaphIsolateDORNs$visit 
StaphIsolateDORNs$subject_timepoint = paste(StaphIsolateDORNs$patient_id, StaphIsolateDORNs$visit,sep="_")
MetaPhlanAbundancetable_Staphylococcus = MetaPhlanAbundancetable_Staphylococcus[c("Staphylococcus_aureus","Timepoint", "Subject", "subject_timepoint")]
Joined = StaphIsolateDORNs %>% left_join(MetaPhlanAbundancetable_Staphylococcus, by="subject_timepoint", drop.)
Joined_noNAs = Joined %>% filter(!is.na(Staphylococcus_aureus))
write.csv(Joined_noNAs, "SAureusAbundances_by_isolateSource.csv")
length(unique(Joined_noNAs$subject_timepoint))
