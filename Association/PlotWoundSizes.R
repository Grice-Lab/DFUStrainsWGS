# Amelia's late Fall 2021 Mouse healing experiments 

WoundSizes = read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/wound_males_females.csv")

WoundSizes = WoundSizes %>% mutate(ConditionType = case_when(group=="PBS" ~ "Control", 
                                                             (group=="LAC" | group=="SA2498") ~ "ReferenceStrain", 
                                                             (group=="ARM72" | group=="ARM81") ~ "Clinical"
                                                           ))


WoundSizes = WoundSizes %>% mutate(ConditionName = case_when(group=="PBS" ~ "Vehicle", 
                                                             group=="LAC" ~ "WildType S. aureus", 
                                                             group=="SA2498" ~ paste0("\u0394", "crt S. aureus"), 
                                                             group=="ARM72" ~ "DORN925 High", 
                                                             group=="ARM81" ~ "DORN1088 Low")
                                                             )

WoundSizes$Sex= if_else(WoundSizes$sex=="m", "Male", "Female")

woundhealplot = ggplot(WoundSizes, aes( group=ConditionName, x=group, y=day_14_change, fill=ConditionType, shape=Sex)) +
  geom_boxplot(alpha=0.5) + geom_jitter( alpha=.9, size=2, width=.1, color="black")
woundhealplot$data$group = factor(woundhealplot$data$ConditionName, levels=c("Vehicle",
                                                                     "DORN925 High",
                                                                     "DORN1088 Low",
                                                                     "WildType S. aureus",
                                                                     paste0("\u0394", "crt S. aureus")))
woundhealplot 


WoundSizesReferences= WoundSizes %>% filter(ConditionType=="ReferenceStrain")
WoundSizesClinical= WoundSizes %>% filter(ConditionType=="Clinical")

wilcRefs = wilcox.test(day_14_change~ConditionName , data=WoundSizesReferences)
wilcRefsP = round(wilcRefs$p.value, digits=5)

wilcClin = wilcox.test(day_14_change~ConditionName , data=WoundSizesClinical)
wilcClinP = round(wilcClin$p.value, digits=5)

wilcControl =  WoundSizes %>% filter(group %in% c("ARM81", "PBS"))
wilc81 = wilcox.test(day_14_change~ ConditionName, data=wilcControl)



woundhealplot+ annotate("text",label=paste0("p=", wilcClinP), x=2.5,y=180) +  annotate("text", label=paste0("p=", wilcRefsP),x=4.5, y=220)+ scale_fill_manual(values=c("darkseagreen", "darkslateblue", "darkorange")) +
  theme_light()  + xlab("S. aureus Strain") + ylab("% Original Wound Size at Day 14")  + ggtitle("Diabetic Mouse Wound Healing") + theme(plot.title=element_text(size=20, hjust=.5))


table(WoundSizes$sex)
