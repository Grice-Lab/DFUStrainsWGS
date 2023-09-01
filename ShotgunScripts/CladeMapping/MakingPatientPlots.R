library(stringr)
library(ggplot2)
library(dplyr)
library(ggpattern)
randpalette18=c("#B300B3","#E6AB02",
                "#0000B3","#006400",
                "#A6761D","#1B9E77",
                "#B3DE69","#FF7F00",
                "#681A1A","#7570B3",
                "#1F78B4","#F2A687",
                "#A6CEE3","#6A3D9A",
                "#666666","#F6BE00",
                "#33A02C","#E6F5C9")

CC_Colorscheme = read.csv("~/Documents/DataInputGithub/data/Phylogeny2022Data/CC_Colorscheme.csv")

# Calculate S. aureus abundance as a % of total *bacteria* in the sample
#####################################################################################
PCtSpecies = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Species_By_Sample_TimeSubj_MetaPhlan.csv")
PCtSpeciesBacterial =  read.csv2("/Users/amycampbell/Documents/DataInputGithub/Ellen/Species_Subj_timepoint_Bacterial_fulltable.tsv", sep="\t")
PCtSpeciesBacterial[,2:(ncol(PCtSpeciesBacterial))] = PCtSpeciesBacterial[,2:(ncol(PCtSpeciesBacterial))] %>% mutate_all(function(x) as.numeric(as.character(x)))

Names=PCtSpeciesBacterial$X

PCtSpeciesBacterial=(PCtSpeciesBacterial %>% select(-X))

for(i in 1:ncol(PCtSpeciesBacterial)){
  PCtSpeciesBacterial[,i] = sapply((PCtSpeciesBacterial[,i]) , function(x) as.numeric(as.character(x)))
  PCtSpeciesBacterial[,i] = 100*(PCtSpeciesBacterial[,i]/sum(PCtSpeciesBacterial[,i]) )
}

BacterialCompositionTransposed = data.frame(t(PCtSpeciesBacterial))
colnames(BacterialCompositionTransposed) = Names


PCtStaph = BacterialCompositionTransposed %>% select(Staphylococcus_aureus)
PCtStaph$SampleID = row.names(PCtStaph)
PCtStaph$Patient_Time = sapply(PCtStaph$SampleID, function(x) (str_split(x,pattern="_sorted_"))[[1]][2])
PCtStaph$Patient_Time = sapply(PCtStaph$Patient_Time, function(x) str_replace(x, pattern="\\.", "_"))
PCtStaph$Patient = sapply(PCtStaph$Patient_Time, function(x) (str_split(x, "_"))[[1]][1])

PCtStaph$Timepoint = sapply(PCtStaph$Patient_Time, function(x) as.numeric((str_split(x, "_"))[[1]][2]) )

PCtStaph$Patient_Time = paste(PCtStaph$Patient, as.character(PCtStaph$Timepoint), sep="_")


# Patient 176
composition176=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_176_composition.csv")
composition176$Patient_Time= paste0("176_", composition176$Timepoint)
composition176 = composition176 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

composition176 %>% select(Timepoint, PctCase)

composition176$CC5 = composition176$PctCase
composition176$CC1 = 100 - composition176$PctCase 

composition176$CC5 = (composition176$CC5*composition176$Staphylococcus_aureus)/100

composition176$CC1 = (composition176$CC1*composition176$Staphylococcus_aureus)/100

MeltedPlot = composition176 %>% select(CC5, CC1, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))
ggplot(composition176, aes(x=Timepoint, y=Staphylococcus_aureus)) + geom_point(color="#B8860B") + geom_line(color="#B8860B") + theme_classic() + scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))

fillColors = (CC_Colorscheme %>% filter(CCLabel %in% c("CC1", "CC5")))$hexval
Patient176Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value, fill=variable ), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus)) + scale_fill_manual(values=fillColors)
Patient176Plot = Patient176Plot + theme_classic() + ylab("% Abundance Staphylococcus aureus") + scale_x_continuous(breaks=sort(composition176$Timepoint)) + labs(fill="Clade")+  scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))
Patient176Plot$data$variable = factor(Patient176Plot$data$variable, levels=c("CC1", "CC5"))
ggsave(Patient176Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient176_UTD.pdf",width=11, height=4)


# Patient 145
composition145=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_145_composition.csv")
composition145$Patient_Time= paste0("145_", composition145$Timepoint)
composition145 = composition145 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

composition145 %>% select(Timepoint, PctCase)

composition145$CC1 = composition145$PctCase
composition145$CC30 = 100 - composition145$PctCase 

composition145$CC30 = (composition145$CC30*composition145$Staphylococcus_aureus)/100

composition145$CC1 = (composition145$CC1*composition145$Staphylococcus_aureus)/100

MeltedPlot = composition145 %>% select(CC30, CC1, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))

fillColors = rev((CC_Colorscheme %>% filter(CCLabel %in% c("CC30", "CC1")))$hexval)

ggplot(composition145, aes(x=Timepoint, y=Staphylococcus_aureus)) + geom_point(color="#B8860B") + geom_line(color="#B8860B") + theme_classic() + scale_x_continuous(breaks=sort(composition145$Timepoint))

Patient145Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value, fill=variable ), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus)) + scale_fill_manual(values=fillColors)
Patient145Plot = Patient145Plot + theme_classic() + ylab("% Abundance Staphylococcus aureus")  + labs(fill="Clade")+  scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))
Patient145Plot$data$variable = factor(Patient145Plot$data$variable, levels=c("CC30", "CC1"))

ggsave(Patient145Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient145_UTD.pdf",width=11, height=4)



# Patient 191
composition191=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_191_composition.csv")
composition191$Patient_Time= paste0("191_", composition191$Timepoint)
composition191 = composition191 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

composition191 %>% select(Timepoint, PctCase)

composition191$CC5 = composition191$PctCase
composition191$CC8 = 100 - composition191$PctCase 

composition191$CC5 = (composition191$CC5*composition191$Staphylococcus_aureus)/100

composition191$CC8 = (composition191$CC8*composition191$Staphylococcus_aureus)/100

fillColors = ((CC_Colorscheme %>% filter(CCLabel %in% c("CC5", "CC8")))$hexval)


MeltedPlot = composition191 %>% select(CC5, CC8, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))
ggplot(composition191, aes(x=Timepoint, y=Staphylococcus_aureus)) + geom_point(color="#B8860B") + geom_line(color="#B8860B") + theme_classic() +  scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))

Patient191Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value, fill=variable ), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus)) + scale_fill_manual(values=fillColors)
Patient191Plot = Patient191Plot + theme_classic() + ylab("% Abundance Staphylococcus aureus") +  scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))+ labs(fill="Clade")
Patient191Plot$data$variable = factor(Patient191Plot$data$variable, levels=c("CC5", "CC8"))
ggsave(Patient191Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient191_UTD.pdf",width=11, height=4)

# Patient 197
composition197=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_197_composition.csv")
composition197$Patient_Time= paste0("197_", composition197$Timepoint)
composition197 = composition197 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

composition197 %>% select(Timepoint, PctCase)

composition197$CC5_1 = composition197$PctCase
composition197$CC5_2 = 100 - composition197$PctCase 

composition197$CC5_1 = (composition197$CC5_1*composition197$Staphylococcus_aureus)/100

composition197$CC5_2 = (composition197$CC5_2*composition197$Staphylococcus_aureus)/100

MeltedPlot = composition197 %>% select(CC5_1, CC5_2, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))
ggplot(composition197, aes(x=Timepoint, y=Staphylococcus_aureus)) + geom_point(color="#B8860B") + geom_line(color="#B8860B") + theme_classic() + scale_x_continuous(breaks=sort(composition197$Timepoint))

# Just temporary; will add two patterns of CC5 color in illustrator
fillColors = (CC_Colorscheme %>% filter(CCLabel%in% c("CC5", "CC9")))$hexval

Patient197Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value,fill=variable), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus))  + scale_fill_manual(values=fillColors)
Patient197Plot = Patient197Plot + theme_classic() + ylab("% Abundance Staphylococcus aureus") + scale_x_continuous(breaks=0:14, limits=c(-.5,14.5)) + labs(fill="Clade")
Patient197Plot$data$variable = factor(Patient197Plot$data$variable, levels=c("CC5_1", "CC5_2"))
ggsave(Patient197Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient197_UTD.pdf",width=11, height=4)


# Patient 173
composition173=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_173_composition.csv")
composition173$Patient_Time= paste0("173_", composition173$Timepoint)
composition173 = composition173 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

composition173 %>% select(Timepoint, PctCase)

composition173$CC72 = composition173$PctCase
composition173$CC5 = 100 - composition173$PctCase 

composition173$CC72 = (composition173$CC72*composition173$Staphylococcus_aureus)/100

composition173$CC5 = (composition173$CC5*composition173$Staphylococcus_aureus)/100
fillColors = rev((CC_Colorscheme %>% filter(CCLabel%in% c("CC5", "CC72")))$hexval)

MeltedPlot = composition173 %>% select(CC72, CC5, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))
ggplot(composition173, aes(x=Timepoint, y=Staphylococcus_aureus)) + geom_point(color="#B8860B") + geom_line(color="#B8860B") + theme_classic() 

Patient173Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value, fill=variable ), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus)) + scale_fill_manual(values=fillColors)
Patient173Plot = Patient173Plot + theme_classic() + ylab("% Abundance Staphylococcus aureus") + scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))  + labs(fill="Clade")
Patient173Plot$data$variable = factor(Patient173Plot$data$variable, levels=c("CC72", "CC5"))
ggsave(Patient173Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient173_UTD.pdf",width=11, height=4)



# Patient 111
composition111=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_111_composition.csv")
composition111$Patient_Time= paste0("111_", composition111$Timepoint)
composition111 = composition111 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

composition111 %>% select(Timepoint, PctCase)

composition111$CC1 = composition111$PctCase
composition111$CC9 = 100 - composition111$PctCase 

composition111$CC1 = (composition111$CC1*composition111$Staphylococcus_aureus)/100

composition111$CC9 = (composition111$CC9*composition111$Staphylococcus_aureus)/100

MeltedPlot = composition111 %>% select(CC1, CC9, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))
ggplot(composition111, aes(x=Timepoint, y=Staphylococcus_aureus)) + geom_point(color="#B8860B") + geom_line(color="#B8860B") + theme_classic() 

fillColors = ((CC_Colorscheme %>% filter(CCLabel%in% c("CC1", "CC9")))$hexval)

Patient111Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value, fill=variable ), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus)) + scale_fill_manual(values=fillColors)
Patient111Plot = Patient111Plot + theme_classic() + ylab("% Abundance Staphylococcus aureus") + scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))  + labs(fill="Clade")
Patient111Plot$data$variable = factor(Patient111Plot$data$variable, levels=c("CC1", "CC9"))
ggsave(Patient111Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient111_UTD.pdf",width=11, height=4)




# Patient 159
composition159=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_159_composition.csv")
composition159$Patient_Time= paste0("159_", composition159$Timepoint)
composition159 = composition159 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

composition159 %>% select(Timepoint, PctCase)

composition159$CC15 = composition159$PctCase
composition159$CC1 = 100 - composition159$PctCase 

composition159$CC15 = (composition159$CC15*composition159$Staphylococcus_aureus)/100

composition159$CC1 = (composition159$CC1*composition159$Staphylococcus_aureus)/100

MeltedPlot = composition159 %>% select(CC15, CC1, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))
ggplot(composition159, aes(x=Timepoint, y=Staphylococcus_aureus)) + geom_point(color="#B8860B") + geom_line(color="#B8860B") + theme_classic() 

fillColors = ((CC_Colorscheme %>% filter(CCLabel%in% c("CC1", "CC15")))$hexval)


Patient159Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value, fill=variable ), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus)) + scale_fill_manual(values=fillColors)
Patient159Plot = Patient159Plot + theme_classic() + ylab("% Abundance Staphylococcus aureus") + scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))  + labs(fill="Clade")
Patient159Plot$data$variable = factor(Patient159Plot$data$variable, levels=c("CC1", "CC15"))
ggsave(Patient159Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient159_UTD.pdf",width=11, height=4)

# Patient 124
composition124=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_124_composition.csv")
composition124$Patient_Time= paste0("124_", composition124$Timepoint)
composition124 = composition124 %>% left_join(PCtStaph %>% select(Patient_Time, Staphylococcus_aureus), by="Patient_Time")

fillColors = rev((CC_Colorscheme %>% filter(CCLabel%in% c("CC5", "CC8", "CC9")))$hexval)

composition124$CC5 = composition124$CC5*composition124$Staphylococcus_aureus/100
composition124$CC8 = composition124$CC8*composition124$Staphylococcus_aureus/100
composition124$CC9 = composition124$CC9*composition124$Staphylococcus_aureus/100

MeltedPlot = composition124 %>% select(CC5, CC9,CC8, Timepoint, Staphylococcus_aureus) %>% reshape2::melt(id.vars=c("Timepoint", "Staphylococcus_aureus"))

Patient124Plot = ggplot(MeltedPlot) + geom_bar(aes(x=Timepoint, y=value, fill=variable ), stat="identity") + geom_line(aes(x=Timepoint, y=Staphylococcus_aureus)) + scale_fill_manual(values=fillColors)

Patient124Plot$data$variable = factor(Patient124Plot$data$variable, levels=rev(c("CC5", "CC8", "CC9")))
Patient124Plot = Patient124Plot  + theme_classic() + ylab("% Abundance Staphylococcus aureus") + scale_x_continuous(breaks=0:14, limits=c(-.5,14.5))  + labs(fill="Clade")
ggsave(Patient124Plot, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/PatientSaureusCompositionPlots/Patient124_UTD.pdf",width=11, height=4)



composition111 = composition111 %>% select(PctCase, Patient_Time)
composition145 = composition145 %>% select(PctCase, Patient_Time)
composition159 = composition159 %>% select(PctCase, Patient_Time)
composition173 = composition173 %>% select(PctCase, Patient_Time)
composition176 = composition176 %>% select(PctCase, Patient_Time)
composition191 = composition191 %>% select(PctCase, Patient_Time)
composition197 = composition197 %>% select(PctCase, Patient_Time)

Total = rbind(composition111,composition145, composition159,composition173,composition176, composition191,composition197  )
TimepointsAll = dim(Total)[1]
Count_All = dim(Total %>% filter(PctCase <10 | PctCase >90))[1]

composition124=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Compositions/Patient_124_composition.csv")
Timepoints124 = dim(composition124)[1]

Count_124 = dim(composition124 %>% filter(CC5 >90 | CC9>90 | CC8 >90))[1]

# 28 / 49
Count_All+Count_124
Timepoints_Total  = TimepointsAll +Timepoints124


CCMap = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Phylogeny2022Data/CCMapPlotting.csv")
PatientInfo = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/DFU_Culture_Collection.csv")
PatientInfo$DORN=paste0("DORN", PatientInfo$Doern.lab.bank.)
PatientInfo = PatientInfo %>% select(DORN, patient_id, visit,)
IncludedIsolates = CCMap %>% filter(!(DORN %in% c("DORN429","DORN1176", "DORN685", "DORN946")))
PatientsCC = IncludedIsolates %>% left_join(PatientInfo, by="DORN")
PatientsCC$Week = PatientsCC$visit*2 
PatientsCC$amount=1
ByCC = PatientsCC %>% group_by(patient_id, CCLabel) %>% summarize(Total=sum(amount)) %>% unique()


listBigTree= c("CC1", "CC12","CC133", "CC15", "CC20", "CC22", "CC30", "CC398", "CC45", "CC5", "CC59", "CC7", "CC72", "CC8", "CC80", "CC9", "CC97" )

Indices_Remove = which(listBigTree %in% setdiff(listBigTree, ByCC$CCLabel))

randpalette14 = randpalette18[-c(3,6,8,15)]

randpalette14
  

WeekHealed$patient_id = WeekHealed$patient
ByCCweek = ByCC %>% left_join(WeekHealed %>% select(patient_id, week_healed)) %>% unique()
ByCCweek = ByCCweek %>% mutate(WeekNumeric = if_else(is.na(week_healed), 100, as.integer(week_healed)))



CCplot = ggplot(ByCCweek, aes(x=factor(patient_id),y=Total, fill=CCLabel)) + geom_bar(stat="identity") + scale_fill_manual(values=randpalette14) + theme_classic() + theme(axis.text.x=element_text(angle=90))


CCplot$data$patient_id = factor(CCplot$data$patient_id , levels= unique((ByCCweek %>% arrange(week_healed))$patient_id))

ggsave(CCplot, width=14,height=4, file="/Users/amycampbell/Documents/Saureus_Genomics_Paper/CCsByPatient.pdf")
# 123 & on did not heal 


