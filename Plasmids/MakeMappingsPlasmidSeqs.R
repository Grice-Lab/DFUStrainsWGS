library(dplyr)

PlasmidDF=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/Filtered_Plasmid_Presence_Absence.csv")
PlasmidSizes=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/FullSummaryPlasmidsSize.csv")
AllPlasmids = colnames(PlasmidDF)[3:ncol(PlasmidDF)]

PathDF = data.frame()
for(plas in 1:length(AllPlasmids)){
  plasmidID = AllPlasmids[plas]
  fastaString=paste0("plasmid_", plasmidID, ".fasta")
  searchID=plasmidID
  if(plasmidID=="NovelPlasmid1"){
    fastaString="plasmid_novel_a631d6c2a3fbb59aed79586f1bc07cf0.fasta"
    searchID="novel_a631d6c2a3fbb59aed79586f1bc07cf0"
  }
  if(plasmidID=="NovelPlasmid2"){
    fastaString="plasmid_novel_fd4748dc670d16378471f93eb5750c70.fasta"
    searchID="novel_fd4748dc670d16378471f93eb5750c70"
  }
  if(plasmidID=="NovelPlasmid3"){
    fastaString="plasmid_novel_7df5c330ea4bd3a4de0d8b8d27ef5f8f.fasta"
    searchID="novel_7df5c330ea4bd3a4de0d8b8d27ef5f8f"
  }
  if(plasmidID=="AE334"){
    fastaString="plasmid_novel_0c4fe61606d2a97e9e987901b0d9d75c.fasta"
    searchID="novel_0c4fe61606d2a97e9e987901b0d9d75c"
  }
  DF = PlasmidSizes %>% filter(primary_cluster_id==searchID)
  IsolateVersion=(DF %>% arrange(-df.size))[1,"IsolateID"]
  pathstring=paste0("/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/", IsolateVersion, "/", fastaString)
  PathDF = rbind(PathDF, c(IsolateVersion, plasmidID,pathstring ))
}
colnames(PathDF) = c("GenomeID", "PlasmidID", "FastaPath")
write.csv(PathDF, file="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/PlasmidFastaPaths.csv")


