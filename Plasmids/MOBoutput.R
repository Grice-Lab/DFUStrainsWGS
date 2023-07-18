# Amy Campbell
# July 2023
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

mobtyperresults=args[1] # 
allGenomeslist=args[2] # path to the full set of included isolates with _Final.fasta extension
OutputFile=args[3]

filelist=list.files(mobtyperresults)
genomelist = list.files(allGenomeslist)

genomenames = sapply(genomelist, function(x) str_split(x,"_")[[1]][1])


FullDF = data.frame()
for(filename in filelist){
  filetoread=(paste0(mobtyperresults,"/",filename))
  df = read.csv2(filetoread, sep="\t")
  genome=str_split(filename,"_")[[1]][1]
  outputDF = data.frame(NearestNeighbor=df$mash_nearest_neighbor,mobility = df$predicted_mobility, primary_cluster_id=df$primary_cluster_id, IsolateID = genome)
  FullDF = rbind(FullDF, outputDF)
}

Presence_Absence_Plasmids = data.frame(matrix(0, length(genomenames), length(unique(FullDF$primary_cluster_id))))
rownames(Presence_Absence_Plasmids) = genomenames
colnames(Presence_Absence_Plasmids) = unique(FullDF$primary_cluster_id)

for(genomeid in unique(FullDF$IsolateID)){
  SmallDF = FullDF %>% filter(IsolateID==genomeid)
  for(i in 1:nrow(SmallDF)){
    clustername=(SmallDF[i,"primary_cluster_id"])
    Presence_Absence_Plasmids[genomeid, clustername] <- 1
  }
}

write.csv(Presence_Absence_Plasmids,file=OutputFile)



