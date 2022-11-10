library(dplyr)
library(stringr)
library(ggplot2)
Blast2GoMappingFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoMappingResults.txt", sep="\t",header=T,row.names=NULL)
Blast2GoAnnotationFile = read.csv2("/Users/amycampbell/b2gWorkspace/Blast2GoAnnotations.txt", sep="\t",header=T,row.names=NULL)
View(Blast2GoAnnotationFile)



dim(Blast2GoMappingFile)
dim(Blast2GoAnnotationFile)


head(Blast2GoAnnotationFile)

length(unique((Blast2GoAnnotationFile %>% filter(!is.na(X.GO)))$SeqName))
Blast2GoAnnotationFileNoComponents = Blast2GoAnnotationFile %>% filter( grepl(GO.IDs, pattern="F:") | grepl(GO.IDs, pattern="P:"))
length(unique((Blast2GoAnnotationFileNoComponents %>% filter(!is.na(X.GO)))$SeqName))




GOIDs_Expanded = Blast2GoAnnotationFile$GO.IDs
length(GOIDs_Expanded)
allGOinstances=c()
for(go in 1:length(GOIDs_Expanded)){
  golist = GOIDs_Expanded[[go]]
  print(length((stringr::str_split(golist, pattern="; "))[[1]]))
  allGOinstances = append(allGOinstances,(stringr::str_split(golist, pattern="; "))[[1]])
}

allGOinstances = allGOinstances[allGOinstances!=""]
GOtable = table(allGOinstances)
GOterm = names(GOtable)
GOdf = data.frame(GOtable)
colnames(GOdf) = c("Term", "Count")
GOdf = GOdf %>% arrange(-Count)
GOdf = GOdf %>% mutate(GO_Type = case_when(grepl(Term, pattern="C:") ~ "Cellular Component", 
                                            grepl(Term, pattern="P:") ~ "Biological Process",
                                            grepl(Term, pattern="F:") ~ "Molecular Function"))
#GOdfNoComponent = GOdf %>% filter(!grepl(Term, pattern="C:"))
GOdfAtLeast5 = GOdf %>% filter(Count >10)
BarPlotGOs = ggplot(GOdfAtLeast5, aes(x=Term, y=Count, fill=GO_Type)) + geom_bar( stat="identity")
BarPlotGOs$data$Term = factor(BarPlotGOs$data$Term, levels=GOdfAtLeast5$Term) 
BarPlotGOs = BarPlotGOs + scale_fill_manual(values=c("#B53737","#190BDA","#E2D629" )) + theme_classic() + theme(axis.text.x=element_text(size=6.5, angle=70, hjust=.4, vjust=.4), plot.title=element_text(size=15, hjust=.5, face="bold")) + ggtitle("GO Terms with >10 Frequency in the Pan-Genome")
#ggsave(BarPlotGOs, width=12,height=6, file="data/Phylogeny2022Data/GOTermsAtLeast10.pdf")


GOdfAtLeast10_Fxns  = GOdfAtLeast5 %>% filter(GO_Type != "Cellular Component")
BarPlotGOsFxns = ggplot(GOdfAtLeast10_Fxns, aes(x=Term, y=Count, fill=GO_Type)) + geom_bar( stat="identity",color="black") + theme_classic() + theme(axis.text.x=element_text(size=8, angle=70, hjust=.4, vjust=.4), plot.title=element_text(size=15, hjust=.5, face="bold")) + scale_fill_manual(values=c("#190BDA","#E2D629" )) + ggtitle("Functional GO Terms with >10 Frequency in the Pan-Genome")
BarPlotGOsFxns$data$Term = factor(BarPlotGOsFxns$data$Term, levels = (GOdfAtLeast10_Fxns %>% arrange(-Count))$Term)

ggsave(BarPlotGOsFxns, width=13,height=7, file="data/Phylogeny2022Data/GOTermsAtLeast10FunctionsOnly.pdf")
