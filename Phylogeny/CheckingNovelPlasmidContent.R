# Amy Campbell
# There were 3 'novel' plasmids (not yet found in MOB-Suite database and also not found in NCBI) 
# that are circularized. However, want to check that their gene contents actually look like plasmids
library(dplyr)
library(stringr)
path_presence_absence = "~/Documents/DataInputGithub/data/IntraPatient/GenePresence_ByPhage_Plasmid.csv"
path_annotationsRoary = "~/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv"
output_plasmid_annotations = "~/Documents/DataInputGithub/data/IntraPatient/Plasmids/AnnotationsNovelPlasmids/"
MOBdatabase = read.csv2("~/Documents/DataInputGithub/data/IntraPatient/Mobsuite_clusters.txt",sep="\t")

presence_absence_by_HGT = read.csv(path_presence_absence)

MOBdatabase = MOBdatabase %>% filter(primary_cluster_id %in% colnames(presence_absence_by_HGT))

savecols = colnames(presence_absence_by_HGT)
GeneAnnotations = read.csv(path_annotationsRoary) %>% select(Gene, Annotation)

presence_absence_by_HGT$FixedGeneNames = sapply(presence_absence_by_HGT$X, function(x) str_replace(x,"\\'", "_"))
presence_absence_by_HGT$FixedGeneNames = sapply(presence_absence_by_HGT$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
presence_absence_by_HGT$FixedGeneNames = sapply(presence_absence_by_HGT$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
presence_absence_by_HGT$FixedGeneNames = sapply(presence_absence_by_HGT$FixedGeneNames, function(x) str_replace(x,":", "_"))
presence_absence_by_HGT$FixedGeneNames = sapply(presence_absence_by_HGT$FixedGeneNames, function(x) str_replace(x,"-", "_"))


GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$Gene, function(x) str_replace(x,"\\'", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,":", "_"))
GeneAnnotations$FixedGeneNames = sapply(GeneAnnotations$FixedGeneNames, function(x) str_replace(x,"-", "_"))


presence_absence_by_HGT = presence_absence_by_HGT %>% left_join(GeneAnnotations,by="FixedGeneNames")


# Novel plasmid 1:
##################
novelplasmid1genes = presence_absence_by_HGT %>% filter(NovelPlasmid1==1) %>% select(Gene, Annotation)
write.csv(novelplasmid1genes, paste0(output_plasmid_annotations, "NovelPlasmid1Genes.csv"))

# Novel plasmid 2:
###################
novelplasmid2genes = presence_absence_by_HGT %>% filter(NovelPlasmid2==1) %>% select(Gene, Annotation)
# contains almost entirely phage-related genes? and also contains sak, scn, etc. 
write.csv(novelplasmid2genes, paste0(output_plasmid_annotations, "NovelPlasmid2_Phage_Genes.csv"))

# Novel plasmid 3:
##################
novelplasmid3genes = presence_absence_by_HGT %>% filter(NovelPlasmid3==1) %>% select(Gene, Annotation)

write.csv(novelplasmid3genes, paste0(output_plasmid_annotations, "NovelPlasmid3_Genes.csv"))

# Checking existing MOB clusters
################################
presence_absence_by_HGT_Updated_Plasmids = presence_absence_by_HGT  %>% select(colnames(presence_absence_by_HGT)[!grepl(x=colnames(presence_absence_by_HGT), pattern="Phage")] )

presence_absence_by_HGT_Updated_Plasmids_phageContainingAnnotations = presence_absence_by_HGT_Updated_Plasmids %>% filter(grepl("phage", Annotation) | grepl("Phage", Annotation))

NumPhageWordGenes_Plasmid = colSums(presence_absence_by_HGT_Updated_Plasmids_phageContainingAnnotations[,2:(ncol(presence_absence_by_HGT_Updated_Plasmids_phageContainingAnnotations)-4)])
potential = NumPhageWordGenes_Plasmid[NumPhageWordGenes_Plasmid>1]


ProbablePhages = c("NovelPlasmid2")

# AF724 
########
AF724Genes = presence_absence_by_HGT %>% filter(AF724==1) #%>% left_join(GeneAnnotations, by="FixedGeneNames")
AF724Genes$Annotation

# the mob cluster contains no relaxase types or rep types
MOBdatabase %>% filter(primary_cluster_id=="AF724")
ProbablePhages = append(ProbablePhages, "AF724")

# AG412
#######
AG412Genes = presence_absence_by_HGT %>% filter(AG412==1)

# the mob cluster contains no relaxase types or rep types
MOBdatabase %>% filter(primary_cluster_id=="AG412")
ProbablePhages = append(ProbablePhages, "AG412")

# AG371
#######
AG371Genes = presence_absence_by_HGT %>% filter(AG371==1)

# the mob cluster contains no relaxase types or rep types
MOBdatabase %>% filter(primary_cluster_id=="AG371")
ProbablePhages = append(ProbablePhages, "AG371")

# AG705
#######
AG705Genes = presence_absence_by_HGT %>% filter(AG705==1)

# the mob cluster contains no relaxase types or rep types
MOBdatabase %>% filter(primary_cluster_id=="AG705")
ProbablePhages = append(ProbablePhages, "AG705")


# AE220
#######
AE220Genes = presence_absence_by_HGT %>% filter(AE220==1)

# the mob cluster contains no relaxase types or rep types
MOBdatabase %>% filter(primary_cluster_id=="AE220")
ProbablePhages = append(ProbablePhages, "AE220")

# AF651
#######
AF651Genes = presence_absence_by_HGT %>% filter(AF651==1)
MOBdatabase %>% filter(primary_cluster_id=="AF651")
ProbablePhages = append(ProbablePhages, "AF651")


presence_absence_by_HGT_Updated = presence_absence_by_HGT %>% select(-ProbablePhages,-FixedGeneNames)
write.csv(presence_absence_by_HGT_Updated, file="~/Documents/DataInputGithub/data/IntraPatient/GenePresence_ByPhage_Plasmid_Updated.csv")




