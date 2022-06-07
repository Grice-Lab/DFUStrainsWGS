
library(tidyr)
library(ggplot2)
library(dplyr)
library(readr)

# Update 6/7/22 
# To include all phenotypes, each subset by KDE on that phenostype (zero-one normalized)

setwd("~/Desktop/GriceLabGit/DFUStrainsWGS/")

# Inflation in Xanthin
#######################
patternsFullXanthin = read_delim("data/dbgwas2022/XanthinZeroOne/bugwasPatternsXanthinZeroOne.txt", " ",col_names=F)
colnames(patternsFullXanthin) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

patternsSubsetXanthin = read_delim("data/dbgwas2022/XanthinSubsetZeroOne/bugwasPatternsXanthinKDEsubset.txt", " ",col_names=F)
colnames(patternsSubsetXanthin) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

pvectFullXanthin =  patternsFullXanthin$`p-value`
pvectSubsetXanthin =  patternsSubsetXanthin$`p-value`


qqobservedFullXanthin= -log10(sort(pvectFullXanthin,decreasing=FALSE))
qqobservedSubsetXanthin = -log10(sort(pvectSubsetXanthin,decreasing=FALSE))

qqexpectFullXanthin = -log10( ppoints(length(pvectFullXanthin) ))
qqexpectSubsetXanthin = -log10( ppoints(length(pvectSubsetXanthin) ))


dataframeFullXanthin= data.frame(expect=qqexpectFullXanthin, observe=qqobservedFullXanthin)
dataframeSubsetXanthin = data.frame(expect=qqexpectSubsetXanthin, observe=qqobservedSubsetXanthin)

plotobjXanthin = ggplot(data=dataframeFullXanthin, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for staphyloxanthin DBGWAS p-values")  + theme_classic() +geom_point(data=dataframeSubsetXanthin, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)

Siderophore 
#############################
patternsFullSiderophore = read_delim("data/dbgwas2022/SiderophoreZeroOne/bugwasPatternsSiderophoreZeroOne.txt", " ",col_names=F)
colnames(patternsFullSiderophore) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

patternsSubsetSiderophore = read_delim("data/dbgwas2022/SiderophoreSubsetZeroOne/bugwasPatternsSiderophoreSubset.txt", " ",col_names=F)
colnames(patternsSubsetSiderophore) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

pvectFullSiderophore =  patternsFullSiderophore$`p-value`
pvectSubsetSiderophore =  patternsSubsetSiderophore$`p-value`


qqobservedFullSiderophore= -log10(sort(pvectFullSiderophore,decreasing=FALSE))
qqobservedSubsetSiderophore = -log10(sort(pvectSubsetSiderophore,decreasing=FALSE))

qqexpectFullSiderophore = -log10( ppoints(length(pvectFullSiderophore) ))
qqexpectSubsetSiderophore = -log10( ppoints(length(pvectSubsetSiderophore) ))


dataframeFullSiderophore= data.frame(expect=qqexpectFullSiderophore, observe=qqobservedFullSiderophore)
dataframeSubsetSiderophore = data.frame(expect=qqexpectSubsetSiderophore, observe=qqobservedSubsetSiderophore)

plotobjsiderophore = ggplot(data=dataframeFullSiderophore, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for Siderophore DBGWAS")  + theme_classic() +geom_point(data=dataframeSubsetSiderophore, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)


# Staphylokinase (sak-positive only)
####################################
patternsFullKinase = read_delim("data/dbgwas2022/staphylokinaseZeroOne/bugwasPatterns_KinaseZeroOne.txt", " ",col_names=F)
colnames(patternsFullKinase) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

patternsSubsetKinase = read_delim("data/dbgwas2022/staphylokinaseSubsetZeroOne/bugwasPatterns_KinaseSubset.txt", " ",col_names=F)
colnames(patternsSubsetKinase) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

pvectFullKinase =  patternsFullKinase$`p-value`
pvectSubsetKinase =  patternsSubsetKinase$`p-value`


qqobservedFullKinase= -log10(sort(pvectFullKinase,decreasing=FALSE))
qqobservedSubsetKinase = -log10(sort(pvectSubsetKinase,decreasing=FALSE))

qqexpectFullKinase = -log10( ppoints(length(pvectFullKinase) ))
qqexpectSubsetKinase = -log10( ppoints(length(pvectSubsetKinase) ))


dataframeFullKinase= data.frame(expect=qqexpectFullKinase, observe=qqobservedFullKinase)
dataframeSubsetKinase = data.frame(expect=qqexpectSubsetKinase, observe=qqobservedSubsetKinase)

plotobjKinase = ggplot(data=dataframeFullKinase, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for Staphylokinase DBGWAS(sak-positive isolates)")  + theme_classic() +geom_point(data=dataframeSubsetKinase, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)
plotobjKinase



# Biofilm
####################################
patternsFullBiofilm = read_delim("data/dbgwas2022/BiofilmZeroOne/bugwasPatternsBiofilmZeroOne.txt", " ",col_names=F)
colnames(patternsFullBiofilm) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

patternsSubsetBiofilm = read_delim("data/dbgwas2022/BiofilmSubsetZeroOne/bugwasPatternsBiofilmSubset.txt", " ",col_names=F)
colnames(patternsSubsetBiofilm) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

pvectFullBiofilm =  patternsFullBiofilm$`p-value`
pvectSubsetBiofilm =  patternsSubsetBiofilm$`p-value`


qqobservedFullBiofilm= -log10(sort(pvectFullBiofilm,decreasing=FALSE))
qqobservedSubsetBiofilm = -log10(sort(pvectSubsetBiofilm,decreasing=FALSE))

qqexpectFullBiofilm = -log10( ppoints(length(pvectFullBiofilm) ))
qqexpectSubsetBiofilm = -log10( ppoints(length(pvectSubsetBiofilm) ))

dataframeFullBiofilm= data.frame(expect=qqexpectFullBiofilm, observe=qqobservedFullBiofilm)
dataframeSubsetBiofilm = data.frame(expect=qqexpectSubsetBiofilm, observe=qqobservedSubsetBiofilm)

plotobjBiofilm = ggplot(data=dataframeFullBiofilm, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for Biofilm DBGWAS")  + theme_classic() +geom_point(data=dataframeSubsetBiofilm, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)



ggsave(gridExtra::grid.arrange(plotobjBiofilm, plotobjXanthin, plotobjKinase, plotobjsiderophore), file="data/dbgwas2022/InflationKDEclusters.pdf")

original = gridExtra::grid.arrange(plotobjBiofilm, plotobjXanthin, plotobjKinase, plotobjsiderophore)


# With the ohter patterns 
##############################
patternsFullXanthin = read_delim("data/dbgwas2022/XanthinZeroOne/step2/patterns.txt", " ",col_names=T)

patternsSubsetXanthin = read_delim("data/dbgwas2022/XanthinSubsetZeroOne/step2/patterns.txt", " ",col_names=T)

pvectFullXanthin =  patternsFullXanthin$`p-value`
pvectSubsetXanthin =  patternsSubsetXanthin$`p-value`


qqobservedFullXanthin= -log10(sort(pvectFullXanthin,decreasing=FALSE))
qqobservedSubsetXanthin = -log10(sort(pvectSubsetXanthin,decreasing=FALSE))

qqexpectFullXanthin = -log10( ppoints(length(pvectFullXanthin) ))
qqexpectSubsetXanthin = -log10( ppoints(length(pvectSubsetXanthin) ))


dataframeFullXanthin= data.frame(expect=qqexpectFullXanthin, observe=qqobservedFullXanthin)
dataframeSubsetXanthin = data.frame(expect=qqexpectSubsetXanthin, observe=qqobservedSubsetXanthin)

plotobjXanthin = ggplot(data=dataframeFullXanthin, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for staphyloxanthin DBGWAS p-values")  + theme_classic() +geom_point(data=dataframeSubsetXanthin, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)

#Siderophore 
#############################
patternsFullSiderophore = read_delim("data/dbgwas2022/SiderophoreZeroOne/step2/patterns.txt", " ",col_names=T)

patternsSubsetSiderophore = read_delim("data/dbgwas2022/SiderophoreSubsetZeroOne/step2/patterns.txt", " ",col_names=T)

pvectFullSiderophore =  patternsFullSiderophore$`p-value`
pvectSubsetSiderophore =  patternsSubsetSiderophore$`p-value`


qqobservedFullSiderophore= -log10(sort(pvectFullSiderophore,decreasing=FALSE))
qqobservedSubsetSiderophore = -log10(sort(pvectSubsetSiderophore,decreasing=FALSE))

qqexpectFullSiderophore = -log10( ppoints(length(pvectFullSiderophore) ))
qqexpectSubsetSiderophore = -log10( ppoints(length(pvectSubsetSiderophore) ))


dataframeFullSiderophore= data.frame(expect=qqexpectFullSiderophore, observe=qqobservedFullSiderophore)
dataframeSubsetSiderophore = data.frame(expect=qqexpectSubsetSiderophore, observe=qqobservedSubsetSiderophore)

plotobjsiderophore = ggplot(data=dataframeFullSiderophore, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for Siderophore DBGWAS")  + theme_classic() +geom_point(data=dataframeSubsetSiderophore, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)


# Staphylokinase (sak-positive only)
####################################
patternsFullKinase = read_delim("data/dbgwas2022/staphylokinaseZeroOne/step2/patterns.txt", " ",col_names=T)

patternsSubsetKinase = read_delim("data/dbgwas2022/staphylokinaseSubsetZeroOne/step2/patterns.txt", " ",col_names=T)

pvectFullKinase =  patternsFullKinase$`p-value`
pvectSubsetKinase =  patternsSubsetKinase$`p-value`


qqobservedFullKinase= -log10(sort(pvectFullKinase,decreasing=FALSE))
qqobservedSubsetKinase = -log10(sort(pvectSubsetKinase,decreasing=FALSE))

qqexpectFullKinase = -log10( ppoints(length(pvectFullKinase) ))
qqexpectSubsetKinase = -log10( ppoints(length(pvectSubsetKinase) ))


dataframeFullKinase= data.frame(expect=qqexpectFullKinase, observe=qqobservedFullKinase)
dataframeSubsetKinase = data.frame(expect=qqexpectSubsetKinase, observe=qqobservedSubsetKinase)

plotobjKinase = ggplot(data=dataframeFullKinase, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for Staphylokinase DBGWAS(sak-positive isolates)")  + theme_classic() +geom_point(data=dataframeSubsetKinase, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)
plotobjKinase



# Biofilm
####################################
patternsFullBiofilm = read_delim("data/dbgwas2022/BiofilmZeroOne/step2/patterns.txt", " ",col_names=T)
#colnames(patternsFullBiofilm) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

patternsSubsetBiofilm = read_delim("data/dbgwas2022/BiofilmSubsetZeroOne/step2/patterns.txt", " ",col_names=T)
colnames(patternsSubsetBiofilm) =c("pattern", "p-value", "q-value", "weight", "wald_statistic")

pvectFullBiofilm =  patternsFullBiofilm$`p-value`
pvectSubsetBiofilm =  patternsSubsetBiofilm$`p-value`


qqobservedFullBiofilm= -log10(sort(pvectFullBiofilm,decreasing=FALSE))
qqobservedSubsetBiofilm = -log10(sort(pvectSubsetBiofilm,decreasing=FALSE))

qqexpectFullBiofilm = -log10( ppoints(length(pvectFullBiofilm) ))
qqexpectSubsetBiofilm = -log10( ppoints(length(pvectSubsetBiofilm) ))

dataframeFullBiofilm= data.frame(expect=qqexpectFullBiofilm, observe=qqobservedFullBiofilm)
dataframeSubsetBiofilm = data.frame(expect=qqexpectSubsetBiofilm, observe=qqobservedSubsetBiofilm)

plotobjBiofilm = ggplot(data=dataframeFullBiofilm, aes(x=expect, y=observe)) + geom_point() +  geom_abline(slope=1, color="red") + xlim(0, 8) + ylim(0,8)  + xlab("-log10(expected p-values)") + ylab("-log10(observed p-values)") +ggtitle("QQ Plot for Biofilm DBGWAS")  + theme_classic() +geom_point(data=dataframeSubsetBiofilm, aes(x=expect, y=observe), color="blue", alpha=.5) #+ geom_point(data=dataframe_uncorrected, aes(x=expect, y=observe), color="black") + geom_abline(slope=1)

new = gridExtra::grid.arrange(plotobjBiofilm, plotobjXanthin, plotobjKinase, plotobjsiderophore)














