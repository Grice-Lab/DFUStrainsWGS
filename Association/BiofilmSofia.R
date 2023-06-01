
library(ggpubr)
DataMay4 = read.csv("Documents/DataInputGithub/data/Validations/SofiaBiofilm/sofia_biofilm_1to1dilution_may4.csv")


DataMay4$Ratio=DataMay4$OD570/DataMay4$OD600

DataMay4Avg = DataMay4 %>% group_by(Strain, Biological) %>% summarize(AverageRatio=mean(Ratio))

MeanJE2= mean((DataMay4Avg %>% filter(Strain=="EGM632"))$AverageRatio)

DataMay4Avg$PctJE2 = 100*(DataMay4Avg$AverageRatio/MeanJE2)
DataMay4Avg = DataMay4Avg %>% filter(Strain!="EGM543")
DataMay4Avg$day="May4"


DataMay12 = read.csv("Documents/DataInputGithub/data/Validations/SofiaBiofilm/sofia_biofilm_1to1dilution_may12.csv")
DataMay12$Ratio=DataMay12$OD570/DataMay12$OD600

DataMay12Avg = DataMay12 %>% group_by(Strain, Biological) %>% summarize(AverageRatio=mean(Ratio))

MeanJE2= mean((DataMay12Avg %>% filter(Strain=="EGM632"))$AverageRatio)

DataMay12Avg$PctJE2 = 100*(DataMay12Avg$AverageRatio/MeanJE2)
DataMay12Avg = DataMay12Avg %>% filter(Strain!="EGM543")
DataMay12Avg$day="May12"


dataBothDays = rbind(DataMay12Avg,DataMay4Avg )


DataMay18 = read.csv("Documents/DataInputGithub/data/Validations/SofiaBiofilm/sofia_biofilm1to1dilution_may18.csv")
DataMay18$Ratio=DataMay18$OD570/DataMay18$OD600

DataMay18Avg = DataMay18 %>% group_by(Strain, Biological) %>% summarize(AverageRatio=mean(Ratio))

MeanJE2= mean((DataMay18Avg %>% filter(Strain=="EGM632"))$AverageRatio)

DataMay18Avg$PctJE2 = 100*(DataMay18Avg$AverageRatio/MeanJE2)
DataMay18Avg$day="May18"


dataAllDays = rbind(dataBothDays,DataMay18Avg )

dataAllDays = dataAllDays%>% filter(Strain!="EGM543")
ggplot(dataAllDays, aes(x=Strain,y=log(PctJE2)))+ geom_boxplot() + geom_jitter(height=0,width=.2,  aes(x=Strain,y=log(PctJE2), color=day)) + stat_compare_means(method="t.test", ref.group="EGM632")

