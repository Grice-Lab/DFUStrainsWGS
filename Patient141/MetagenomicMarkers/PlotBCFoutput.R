# Amy Campbell
# Plot read depth at markers under different conditions
# False positives 
# Low xanthin true positives
# DORN1000, DORN1088, DORN1038


# High xanthin true positives 
# DORN882
library(dplyr)
library(ggplot2)

ParameterRange = read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/ParameterRangeUTD.csv")
ParameterRangeN0 = read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/ParameterRange_N0s.csv")

unique(ParameterRange$Source)

# For N0, for each set of parameters, we have:

#     High true positives: DORN1061, DORN925
#     Low true positives : DORN1088, DORN1000 --> together these make up four groups contributing to true positives

#     False positives: StaphEpi, StaphPettenkoferi

# for each set of parameters and each marker (AllParams), make two dataframes, summarizing each group's depth for each AllParam as:
#           AverageDepth = mean(totalUseableDepth)
#           MinDepth = min(totalUseableDepth)
#           MaxDepth = max(totalUseableDepth)


ParameterRangeN0$AllParams=paste(ParameterRangeN0$Dvalue,ParameterRangeN0$Rvalue, ParameterRangeN0$Lvalue, ParameterRangeN0$Qvalue,ParameterRangeN0$MarkerID, ParameterRangeN0$S2value, ParameterRangeN0$Nvalue,  sep="_")
ParameterRangeN0 = ParameterRangeN0 %>% mutate(Group=case_when( (Source %in% c("DORN1061", "DORN925", "DORN1088","DORN1000")) ~ "TruePositive",
                                                            (Source %in% c("StaphEpi", "StaphPettenkoferi")) ~ "FalsePositive", 
                                                            TRUE ~ "RealData"))

ParameterRangeN0False = ParameterRangeN0 %>% filter(Group=="FalsePositive")
ParameterRangeN0False = ParameterRangeN0False %>% arrange(AllParams)
ParameterRangeN0FalseJustParams = ParameterRangeN0False %>% select(AllParams, MarkerID) %>% unique()

FalseMeans = ParameterRangeN0False %>% group_by(AllParams) %>% summarize(MeanDepth = mean(TotalUseableDepth))
FalseMins = ParameterRangeN0False %>% group_by(AllParams) %>% summarize(MinDepth = min(TotalUseableDepth))
FalseMaxs = ParameterRangeN0False %>% group_by(AllParams) %>% summarize(MaxDepth = max(TotalUseableDepth))

N0sFalse = data.frame(Params=ParameterRangeN0FalseJustParams$AllParams, Marker=ParameterRangeN0False$MarkerID, FalseMean = FalseMeans$MeanDepth, FalseMin=FalseMins$MinDepth, FalseMax= FalseMaxs$MaxDepth)

ParameterRangeN0True = ParameterRangeN0 %>% filter(Group=="TruePositive")
ParameterRangeN0True = ParameterRangeN0True %>% arrange(AllParams)
ParameterRangeN0TrueJustParamsMarkers = ParameterRangeN0True %>% select(AllParams, MarkerID) %>% unique()

TrueMeans = ParameterRangeN0True %>% group_by(AllParams) %>% summarize(MeanDepth = mean(TotalUseableDepth))
TrueMins = ParameterRangeN0True %>% group_by(AllParams) %>% summarize(MinDepth = min(TotalUseableDepth))
TrueMaxs = ParameterRangeN0True %>% group_by(AllParams) %>% summarize(MaxDepth = max(TotalUseableDepth))

N0sTrue = data.frame(Params=ParameterRangeN0True$AllParams,Marker=ParameterRangeN0True$MarkerID, TrueMean=TrueMeans$MeanDepth, TrueMax=TrueMaxs$MaxDepth, TrueMin=TrueMins$MinDepth)

N0sDF = N0sFalse %>% left_join(N0sTrue, by="Params")

ggplot(data = N0sDF,aes(x = FalseMean ,y = TrueMean))  + geom_point() + ggtitle("Depth of False Positive vs. True Positive Reads Matched by BT Parameters -- N -0")#+ geom_errorbar(aes(ymin=TrueMin, ymax=TrueMax),color="lightgray") + geom_errorbarh(aes(xmin=FalseMin, xmax=FalseMax), color="lightgray") + theme_minimal()

result =N0sDF %>% filter(TrueMean==max(N0sDF$TrueMean) & FalseMean == min(N0sDF$FalseMean))

resultOriginal=ParameterRangeN0True %>% filter(AllParams %in% result$Params)

resultOriginal1 = resultOriginal %>% filter(MarkerID=="Marker1")
(resultOriginal1 %>% arrange(-Dvalue, -Rvalue, Lvalue, S2value))[1,]


# Filtered 141_04: filtered_sorted_141-04_30_9_14_05_N0_0.bcf
##################
# Marker 1 depth DP=15;I16=0,0,8,7 # same in N1 version

# Filtered 141_01: filtered_sorted_141-01_30_9_14_05_N0_0.bcf
# Marker 1 Depth DP=11;I16=3,8,0,0 # same in N1 version




resultOriginal2 = resultOriginal %>% filter(MarkerID=="Marker2")
(resultOriginal2 %>% arrange(-Dvalue, -Rvalue, Lvalue, S2value))[1,]
# 30_9_14_10_N0_1
# filtered_sorted_141-04_30_9_14_10_N0_1.bcf
# Marker 2 depth = 0 (uncovered) # Same in N1 version

# filtered_sorted_141-01_30_9_14_10_N0_1.bcf
# Marker 2 depth DP=1;I16=0,1,0,0 # same in N1 version

# S. epi and S. petten still 0 in N1 version 

resultOriginal3 = resultOriginal %>% filter(MarkerID=="Marker3")
(resultOriginal3 %>% arrange(-Dvalue, -Rvalue, Lvalue, S2value))[1,]
#  30_9_18_05_N0_5     

# filtered_sorted_141-01_30_9_18_05_N0_5.bcf
# Marker 3 depth DP=6;I16=4,2,0,0  # same in N1 version

# filtered_sorted_141-04_30_9_18_05_N0_5.bc
# Marker 3 depth DP=1;I16=0,0,0,1 # same in N1 version

# BUT S. epi has DP=9;I16=4,5, 0, 0 in N1 version of Marker 3 (S. petten still 0)


# Total Depth filtered_sorted_141-01
# High xanthin allele:
#  11 + 1 + 6 
# Low xanthin allele: 
#  0 + 0 + 0 


# Total Depth filtered_sorted_141-04
# 
# Low xanthin allele: 
# 15 + 0 + 1
# High Xanthin allele:
# 0 + 0 + 0 


# Filtered_sorted_141_01
# 17/17 high xanthin 

# Filtered_sorted_141_04 
# 16/16 low xanthin


# Marker 1: 
# filtered_sorted_<ID>_30_9_14_05_N0_0.bcf

# Marker 2:
# filtered_sorted_141-04_30_9_14_10_N0_1.bcf

# Marker 3:
# filtered_sorted_<ID>_30_9_18_05_N0_5.bcf



# Trying with N1:



ParameterRange1000 = ParameterRange %>% filter(Source=="DORN1000") %>% select(Source, AllParams, TotalUseableDepth, MarkerID)
ParameterRange1000$Type = "Positive_Low"

ParameterRange1061 = ParameterRange %>% filter(Source=="DORN1061") %>% select(Source, AllParams, TotalUseableDepth, MarkerID)
ParameterRange1061$Type = "Positive_High"

ParameterRange1088 = ParameterRange %>% filter(Source=="DORN1088") %>% select(Source, AllParams, TotalUseableDepth, MarkerID)
ParameterRange1088$Type= "Positive_Low"

ParameterRange925 = ParameterRange %>% filter(Source=="DORN925") %>% select(Source, AllParams, TotalUseableDepth, MarkerID)
ParameterRange925$Type = "Postiive_High"

ParameterRangeSEpi = ParameterRange %>% filter(Source=="StaphEpi") %>% select(Source, AllParams, TotalUseableDepth)
ParameterRangeSEpi$Type = "Negative_Epi"

ParameterRangeSPett = ParameterRange %>% filter(Source=="StaphPettenkoferi") %>% select(Source, AllParams, TotalUseableDepth)
ParameterRangeSPett$Type = "Negative_Petten"

# Iterate through the different positive cases
##############################################

# DORN1000
Compare1000_SPett = ParameterRange1000 %>% left_join(ParameterRangeSPett, by="AllParams")
colnames(Compare1000_SPett) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")
Compare1000_SEpi = ParameterRange1000 %>% left_join(ParameterRangeSEpi, by="AllParams")
colnames(Compare1000_SEpi) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")

# DORN1061
Compare1061_SPett = ParameterRange1061 %>% left_join(ParameterRangeSPett, by="AllParams")
colnames(Compare1061_SPett) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")
Compare1061_SEpi = ParameterRange1061 %>% left_join(ParameterRangeSEpi, by="AllParams")
colnames(Compare1061_SEpi) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")

# DORN925
Compare925_SPett = ParameterRange925 %>% left_join(ParameterRangeSPett, by="AllParams")
colnames(Compare925_SPett) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")
Compare925_SEpi = ParameterRange925 %>% left_join(ParameterRangeSEpi, by="AllParams")
colnames(Compare925_SEpi) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")

# DORN1088
Compare1088_SPett = ParameterRange1088 %>% left_join(ParameterRangeSPett, by="AllParams")
colnames(Compare1088_SPett) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")
Compare1088_SEpi = ParameterRange1088 %>% left_join(ParameterRangeSEpi, by="AllParams")
colnames(Compare1088_SEpi) = c("SourceTP", "Parameters", "TotalDepthTP","Marker", "TypeTP", "SourceFP", "DepthFP", "TypeFP")

# 1000, 1061, 1088, 925
FourColorPalette=c("#FFCB3E", "#C1549C","#FB836F", "#7E549F")
ggplot(Compare1000_SPett) + geom_point(color="#FFCB3E", aes(x=TotalDepthTP, y=DepthFP, shape=sourceFP), shape=19, alpha=.5)+
  geom_point(data=Compare1000_SEpi, color="#FFCB3E", aes(x=TotalDepthTP, y=DepthFP), shape=17,alpha=.5) + 
  geom_point(data=Compare1061_SPett, color="#C1549C", aes(x=TotalDepthTP, y=DepthFP), shape=19,alpha=.5) +
  geom_point(data=Compare1061_SEpi, color="#C1549C", aes(x=TotalDepthTP, y=DepthFP), shape=17,alpha=.5) +
  geom_point(data=Compare1088_SPett, color="#FB836F", aes(x=TotalDepthTP, y=DepthFP), shape=19,alpha=.5) +
  geom_point(data=Compare1088_SEpi, color="#FB836F", aes(x=TotalDepthTP, y=DepthFP), shape=17,alpha=.5) +
  geom_point(data=Compare925_SPett, color="#7E549F", aes(x=TotalDepthTP, y=DepthFP), shape=19,alpha=.5) +
  geom_point(data=Compare925_SEpi, color="#7E549F", aes(x=TotalDepthTP, y=DepthFP), shape=17,alpha=.5) + facet_grid(.~Marker)





Depths_PositiveHigh = rbind(ParameterRange925,ParameterRange1061 ) %>% arrange(AllParams)
Depths_PositiveLow= rbind(ParameterRange1000,ParameterRange1088 ) %>% arrange(AllParams)

Depths_PositiveBoth = rbind(Depths_PositiveLow, Depths_PositiveHigh) %>% arrange(AllParams)

Depths_negative = rbind(ParameterRangeSEpi,ParameterRangeSPett ) %>% arrange(AllParams)
Depths_negative2 = rbind(Depths_negative, Depths_negative) %>% arrange(AllParams)


positives_negatives=data.frame(TruePositives=Depths_PositiveBoth$TotalUseableDepth, FalsePositives=Depths_negative2$TotalUseableDepth, negativesParams=Depths_negative2$AllParams, positivesParams=Depths_PositiveBoth$AllParams, sourceFP=Depths_negative2$Source, sourceTP=Depths_PositiveBoth$Source)

ggplot(positives_negatives, aes(x=TruePositives, y=FalsePositives, color=sourceTP, shape=sourceFP)) + geom_point() + scale_color_manual(values=FourColorPalette)


positives_negatives
View(positives_negatives %>% arrange(FalsePositives, -TruePositives))
