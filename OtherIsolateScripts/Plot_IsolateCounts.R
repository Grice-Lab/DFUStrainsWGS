# Amy Campbell
# Generate visualization of which patients 
# at which timepoints we have isolates for 

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
isolates = read.csv("/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/DFU_Staph_aureus_isolates.csv")
patient_metadata = read.csv("/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/gardner_metadata23DEC1.csv")
sixteen_s = read.csv("/Users/amycampbell/Desktop/Club_Grice/scripts/loesche/dfu_uclust/metadata/v1v3_sample_key.csv")
loesche_data = read.table("/Users/amycampbell/Desktop/Club_Grice/scripts/acampbe/DFU/HMMUFOtu_Output/MiSeqV1V3_22_hmmufotu_OTUTable.txt", sep = "\t", header=TRUE)
paperdata = read.csv("/Users/amycampbell/Desktop/Club_Grice/metamap.csv")
paperdata$SampleID = paste(paperdata$patient_id, paperdata$visit, sep="-")


present=colnames(loesche_data)
#setdiff(patient_metadata$patient_id, unique(sixteen_s$patient))
present = sapply(present, function(x) sub("X","", x))

S = sixteen_s %>% filter(SampleID %in% present)
length(unique(patient_metadata$studyid))

#length(unique(paper_data$patient_id))

sixteen = subset(S, SubjectID!="Control")
sixteen$patient = sapply(sixteen$SubjectID, function(x) substr(x, 1, 3))

Presumed_Included = (intersect(patient_metadata$studyid, sixteen$patient))

patient_metadata$patient_id = patient_metadata$studyid

Sfinal = sixteen %>% filter(patient %in% Presumed_Included)

patient_metadata_subset = patient_metadata %>% filter(studyid %in% Presumed_Included)

###########################################################################################
# NOTE : Previous version of this script still categorized 117 and 312 as unhealed. To 
#        restore this, you have to change "study_week_healed" back to wks_to_heal
###########################################################################################
patient_metadata_subset = patient_metadata %>% mutate(healed=ifelse(visit_healed > 6 | is.na(study_week_healed), "Unhealed", "Healed"))

Timepoints = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")
isolates$SampleID = paste((isolates$patient_id), (isolates$visit), sep="-")


length(Timepoints)

subj_timepoint_vector = c()
for (t in Timepoints) {
  for (s in Presumed_Included) {
    subj_timepoint_vector=c(subj_timepoint_vector, paste((s),(t), sep="-"))
  }
}
newdf = data.frame(subj_timepoint_vector)
colnames(newdf) = c("SampleID")


# a simplified isolates frame that has "SampleID", "patient_id", "visit", "IsolateCount" where isolatecount is how many isolates we have from that sample (subject/timepoint combo)
samples_metagenomics = paperdata$SampleID

isolates_simplified = isolates %>% group_by(SampleID) %>% tally() 
isolates_simplified = data.frame(isolates_simplified)
isolates_simplified = newdf %>% left_join(isolates_simplified, by="SampleID")
isolates_simplified = isolates_simplified %>% separate(SampleID, sep="-", into=c("subject", "visit"), remove=FALSE)
isolates_simplified = data.frame(isolates_simplified)
isolates_simplified[is.na(isolates_simplified)] <- 0
isolates_simplified = isolates_simplified %>% mutate(ShotgunPresent = ifelse(SampleID %in% samples_metagenomics, "TRUE", "FALSE"))
isolates_simplified = isolates_simplified %>% group_by(subject) %>% mutate(anyshotgun = ifelse(ShotgunPresent==TRUE, 1, 0))
isolates_simplified = isolates_simplified %>% group_by(subject) %>% filter(sum(n) > 0 | sum(anyshotgun) > 0)
isolates_simplified$visit = sapply(isolates_simplified$visit, function(x) as.integer(x))

isolates_simplified = data.frame(isolates_simplified)
#isolates_simplified$subject = as.factor(isolates_simplified$subject)

# In timepoints 0 through 2 with >0 isolates and metagenomic shotgun present
double_whammy = isolates_simplified %>% filter(ShotgunPresent==TRUE & visit <=2 & n>0)
length(unique(double_whammy$subject))
write.csv(double_whammy, "Looking_For_AGR-C.csv")
subjects_metagenomecounts = data.frame(isolates_simplified %>% group_by(subject) %>% summarise(Metagenomes=sum(ShotgunPresent==TRUE)))
subjects_isolate_timepoint_counts = data.frame(isolates_simplified %>% group_by(subject) %>% summarise(IsolateTimepoints = sum(n>0)))

Counts_Isolates = subjects_metagenomecounts %>% full_join(subjects_isolate_timepoint_counts, by="subject")
with(new, table(ShotgunPresent, subject))
#factor(isolates_simplified$subject, levels = isolates_simplified$subject[order(isolates_simplified$hasany)])

Just_Isolates = Counts_Isolates %>% filter(IsolateTimepoints>3 & Metagenomes==0) 
Just_Metagenomes = Counts_Isolates  %>% filter(IsolateTimepoints==0 & Metagenomes>3) 

sort(table((isolates_simplified %>% filter(n>0))$subject))
subjects_with_atleast4isolateTimepoints <- c(124,149,171, 178, 183, 106, 159, 173, 197, 201, 191, 132, 141, 166, 176, 186 )
#order(isolates_simplified, hasany)
isolate_count_palette = c("#FFFFFF", "lightgoldenrod1",  rev(RColorBrewer::brewer.pal(9, name="RdYlBu")[c(1,2,3)]))
  #c("#FFFFFF", rev(RColorBrewer::brewer.pal(9, name="Spectral")[c(1, 2, 4, 5)])) #RColorBrewer::brewer.pal(9, name="YlGn")[c(1, 4, 6, 8, 9)]

isolateplot <- ggplot(isolates_simplified, aes(x=subject, y=visit, fill = factor(n))) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=12),panel.background=element_rect(fill="lightgrey"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5, size=20), axis.text.x=element_text(angle=270, size=15), axis.text.y=element_text(size=15), axis.title = element_text(size=15)) +
  geom_tile(aes(fill = factor(n), color=factor(ShotgunPresent), width=.9, height=.9), size=.7) +
  scale_fill_manual(values=isolate_count_palette)  + labs(fill='# S. aureus Isolates', color='Metagenomic Shotgun \n Profile Present', title="S. aureus Isolates & Metagenomic Data Collected by Subject/Timepoint") +
  scale_colour_manual(values=c("white", "black")) + scale_y_continuous(breaks=0:13, labels=Timepoints) +
  guides(color = guide_legend(override.aes = list(size=2, stroke=.5, fill=NA)))  + xlab("Subject ID") + ylab("Visit #") # + ggtitle("S. aureus Isolates & Metagenomic Data Present by Subject/Timepoint")
#isolateplot
postscript("IsolatesByTimeSubject.eps", height=10, width=50)
print(isolateplot)
dev.off()



#patient_metadata_subset = patient_metadata %>% filter(studyid %in% isolates_simplified$subject)
#patient_metadata_subset = patient_metadata %>% mutate(healed=ifelse(wks_to_heal > 12.0 | is.na(wks_to_heal), "Unhealed", "Healed"))
healingstatus = patient_metadata_subset[c("healed", "patient_id")]
colnames(healingstatus) = c("healed", "subject")
healingstatus$subject = sapply(healingstatus$subject, function(x) toString(x))
isolates_simplified1 = healingstatus %>% inner_join(isolates_simplified, by="subject")
dim(isolates_simplified1)
#+ geom_tile(data=isolates_simplified,aes(x=subject, y=visit, fill=NA,),size=2,fill=NA)
#theme(legend.title= element_text("Isolate Present"))
length(unique(isolates_simplified$subject))
healeddf = isolates_simplified1 %>% filter(healed =="Healed")
unhealeddf = isolates_simplified1 %>% filter(healed =="Unhealed")
#length(unique(healeddf$subject)) +length(unique(unhealeddf$subject)) 
levels(as.factor(healeddf$n)) 
healeddf$n <- factor(healeddf$n)
levels(healeddf$n) = c("0", "1", "2", "3", "4")
  
unhealed_subset = unhealeddf %>%arrange(n) %>% filter(ShotgunPresent==TRUE)
zeros_unhealed = setdiff(unique(unhealeddf$subject), unique(unhealed_subset$subject))
unhealed_subset= unhealed_subset %>% group_by(subject) %>% summarise(last = max(visit)) %>% arrange(last)

desired_order = append(zeros_unhealed, unhealed_subset$subject)
unhealedplot <- ggplot(unhealeddf, aes(x=subject, y=visit, fill = factor(n))) +
  theme(legend.title =element_text(size=18),legend.text=element_text(size=12), panel.background=element_rect(fill="lightgrey"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5, size=20), axis.text.x=element_text(angle=270, size=15), axis.text.y=element_blank(), axis.ticks.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_tile(aes(fill = factor(n), color=factor(ShotgunPresent), width=.9, height=.9), size=.7) +
  scale_fill_manual(values=isolate_count_palette)  + labs(fill='# S. aureus Isolates', color='Metagenomic Shotgun \n Profile Present', title="Unhealed Samples") +
  scale_colour_manual(values=c("white", "black"))+ xlab("") + scale_y_continuous(breaks=0:13, labels=Timepoints) +
  guides(color = guide_legend(override.aes = list(size=5, stroke=.5, fill=NA))) +  theme(plot.title = element_text(vjust = -8, size=12), axis.text.x=element_text(angle=270, size=12, vjust=.5), axis.title.x = element_blank())   # + ggtitle("S. aureus Isolates & Metagenomic Data Present by Subject/Timepoint")
unhealedplot$data$subject = factor(unhealedplot$data$subject , levels =desired_order)


unhealedplot
# hacky workaround -- theres no n=3 in the levels of healed samples counts, so I'm getting rid of the current '3' color in favor of the '4' color, which it does have 
healed_subset = healeddf %>%arrange(n) %>% filter(ShotgunPresent==TRUE) 
zeros_healed = setdiff(unique(healeddf$subject), unique(healed_subset$subject))
healed_subset = healed_subset %>% group_by(subject) %>% summarise(last = max(visit)) %>% arrange(last)
desired_order = append(zeros_healed, healed_subset$subject)
palette_reduced = c( "#FFFFFF","lightgoldenrod1", "#FDAE61", "#D73027")
healedplot <- ggplot(healeddf, aes(x=subject, y=visit, fill = factor(n))) +
  theme(legend.position = "none", panel.background=element_rect(fill="lightgrey"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5, size=20), axis.text.x=element_text(angle=270, size=15), axis.text.y=element_text(size=15), axis.title = element_text(size=15)) +
  geom_tile(aes(fill = factor(n), color=factor(ShotgunPresent), width=.9, height=.9), size=.7) +
  scale_fill_manual(values=palette_reduced)  + 
  scale_colour_manual(values=c("white", "black")) + labs(title="Healed Samples")  + ylab("Visit #") + scale_y_continuous(breaks=0:13, labels=Timepoints) + 
  theme(plot.title = element_text(vjust = -8, size=12), axis.text.x=element_text(angle=270, size=12, vjust=.5), axis.title.x = element_blank()) 
  #guides(color = guide_legend(override.aes = list(size=2, stroke=.5, fill=NA)))  + xlab("Subject ID") + ylab("Visit #") # + ggtitle("S. aureus Isolates & Metagenomic Data Present by Subject/Timepoint")
healedplot$data$subject = factor(healedplot$data$subject , levels =desired_order)


grid.arrange(healedplot,unhealedplot, nrow=1, ncol=2, widths=c(.8, .6), top=textGrob("S. aureus Isolates & Metagenomic Data Collected by Subject/Timepoint", gp=gpar(fontsize=20, fontface="bold")), bottom=textGrob("Subject ID"))

