
library(tidyr)
library(ggplot2)
library(dplyr)
library(readr)

# data is in /Users/amycampbell/Documents/TestingXanthinDBGWASSubsets
# where, on the LPC:
# Subset219Unitigs.tsv = /home/acampbe/DFU/data/WGS_2020/DBGWAS_Test_219_Xanthin/textualOutput/all_comps_nodes_info.tsv
# Subset104Unitigs.tsv = /home/acampbe/DFU/data/WGS_2020/DBGWAS_Test_104_Xanthin/textualOutput/all_comps_nodes_info.tsv
# Subset98Unitigs.tsv = /home/acampbe/DFU/data/WGS_2020/DBGWAS_Test_98_Xanthin/textualOutput/all_comps_nodes_info.tsv
# BUGWASPatterns219.txt = /home/acampbe/DFU/data/WGS_2020/DBGWAS_Test_219_Xanthin/step2/patterns.txt
# BUGWASPatterns104.txt = /home/acampbe/DFU/data/WGS_2020/DBGWAS_Test_104_Xanthin/step2/patterns.txt
# BUGWASPatterns98.txt = /home/acampbe/DFU/data/WGS_2020/DBGWAS_Test_98_Xanthin/step2/patterns.txt


setwd("/Users/amycampbell/Documents/TestingXanthinDBGWASSubsets")


# Inflation in each
###################

# Full set of isolates Amelia's using
bugWASpatterns219 = read_delim("BUGWASPatterns219.txt", " ", col_names=T)

# K means clusters + patients + 50 SNP cutoff 
bugWASpatterns104 = read_delim("BUGWASPatterns104.txt", " ", col_names=T)

# K means clusters + patients + MLST CCs 
bugWASpatterns98 =read_delim("BUGWASPatterns98.txt", " ", col_names=T)

pvect219 =  bugWASpatterns219$`p-value`
pvect104 =  bugWASpatterns104$`p-value`
pvect98 = bugWASpatterns98$`p-value`

qqobserved219 = -log10(sort(pvect219,decreasing=FALSE))
qqobserved104 = -log10(sort(pvect104,decreasing=FALSE))
qqobserved98 = -log10(sort(pvect98,decreasing=FALSE))

qqexpect219 = -log10( ppoints(length(pvect219) ))
qqexpect104 = -log10( ppoints(length(pvect104) ))
qqexpect98 = -log10( ppoints(length(pvect98) ))


dataframe219 = data.frame(expect=qqexpect219, observe=qqobserved219)
dataframe104 = data.frame(expect=qqexpect104, observe=qqobserved104)
dataframe98 = data.frame(expect=qqexpect98, observe=qqobserved98)

plotobj = ggplot(data=dataframe219, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for DBGWAS p-values")  + theme_classic()#+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)
plotobj= plotobj + geom_point(data=dataframe104, aes(x=expect, y=observe), color="blue", alpha=.5) + geom_point(data=dataframe98, aes(x=expect, y=observe), color="green", alpha=.5) +
annotate(geom="text",color="green", label="98-isolate subset", x=1, y=7, size=8) + annotate(geom="text",color="blue", label="104-isolate subset", x=1, y=7.5, size=8) + annotate(geom="text",color="black", label="219-isolate set", x=1, y=6.5, size=8)



results219 = read_tsv("Subset219Unitigs.tsv")
results104 = read_tsv("Subset104Unitigs.tsv")

results219_top20 = results219 %>% select(`p-value`, `q-Value`,Sequence) %>% arrange(`q-Value`,`p-value` )
results219_top20 = results219_top20[1:20, ]

# Selected top 47 bc there are two top q-values, when you sort by q then p to include all the ones sharing the p (so probably equally distributed)
results104_top20 = results104 %>% select(`p-value`, `q-Value`,Sequence) %>% arrange(`q-Value`, `p-value`)
results104_top44 = results104_top20[1:44, ]

IntersectingUnitigs = intersect(results219_top20$Sequence,results104_top47$Sequence)
View(results104 %>% filter(Sequence %in% IntersectingUnitigs))

ggsave(plotobj, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Plots_Output/InflationPlotThreeClusteringTypes.png")

# a7185595174bd26cbba955f405c06990_1798_arcC
# "bc50a8e1ddf06f071ff651359644e735_2682_group_1096~~~a7185595174bd26cbba955f405c06990_1344_extdb:pgaptmp_001653"  "Select seq ref|WP_000781347.1|	superantigen-like protein SSL5 [Staphylococcus aureus]"
# "a7185595174bd26cbba955f405c06990_1937_extdb:pgaptmp_002241~~~00f10a58fa7304043b40a4fcde135512_2886_extdb:pgaptmp_002971" "Select seq tpg|HBC4173420.1|	TPA: gluconate permease [Staphylococcus aureus] Select seq ref|WP_031897099.1|	gluconate:H+ symporter [Staphylococcus aureus]  which I htink is gntP https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1797221/
# "a7185595174bd26cbba955f405c06990_270_extdb:pgaptmp_000486" shikimate kinase [Staphylococcus aureus]
# And two NAs

# 

#NAs:
 # TATTTTCACAAAATACTATAATGAGGATAGTAAATAGAGAGGAG
