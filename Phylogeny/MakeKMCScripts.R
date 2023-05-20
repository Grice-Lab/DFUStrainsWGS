# Amy Campbell
# Making operations files and shell scripts for each CC I want unique & ubiquitous Kmers for
library(dplyr)
library(stringr)


CCReps = read.csv("Documents/DataInputGithub/CladeRepresentatives.csv")
OutputFolder = "Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/KMCscripts/"
OutputKmers = "/home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/"
KmersPrefixLPC = "/home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/Kmers/"
CCReps = CCReps %>% mutate(KmerName=if_else(grepl("DORN", GenomeID), paste0(GenomeID, "_Final_Kmers"), paste0(GenomeID, "_Kmers")))


for(cc in unique(CCReps$CC)){
  IngroupRows = CCReps %>% filter(CC==cc)
  OutgroupRows = CCReps %>% filter(CC!=cc)
  
  # Operations file for intersecting all the 'in-group' genomes' kmers 
  IngroupIntersectOps = paste0(OutputFolder, cc, "_ingroups.txt")
  
  LineListIntersect = c("INPUT:")
  IngroupIDlist = c()
  for(rowid in 1:nrow(IngroupRows)){
    GroupLine = paste0("Ingroup", rowid, " = ", KmersPrefixLPC,CCReps[rowid, 4] )
    LineListIntersect=append(LineListIntersect, GroupLine)
    IngroupIDlist=append(IngroupIDlist, paste0("Ingroup", rowid))
  }
  LineListIntersect = append(LineListIntersect, "")
  LineListIntersect = append(LineListIntersect, "OUTPUT:")
  CommandLine = paste0(OutputKmers, cc, "_intersect = ")
  CommandLine=paste0(CommandLine, str_c(IngroupIDlist, collapse = " * "))
  
  LineListIntersect = append(LineListIntersect, CommandLine)
  
  writeLines(LineListIntersect, IngroupIntersectOps)
  
  
  
  # Operations file for finding union between all the 'out-group' genomes' kmers 
  OutgroupUnionOps = paste0(OutputFolder, cc, "_outgroups.txt")
  LineListUnion = c("INPUT:")
  OutgroupIDlist = c()
  for(rowid in 1:nrow(OutgroupRows)){
    GroupLine = paste0("Outgroup", rowid, " = ", KmersPrefixLPC,CCReps[rowid, 4] )
    LineListUnion=append(LineListUnion, GroupLine)
    OutgroupIDlist=append(OutgroupIDlist, paste0("Outgroup", rowid))
  }
  LineListUnion = append(LineListUnion, "")
  LineListUnion = append(LineListUnion, "OUTPUT:")
  CommandLine = paste0(OutputKmers, cc, "_Union = ")
  CommandLine=paste0(CommandLine, str_c(OutgroupIDlist, collapse = " + "))
  
  LineListUnion = append(LineListUnion, CommandLine)
  
  writeLines(LineListUnion, OutgroupUnionOps)
  
  
  # Shell script
  IntersectOutput=paste0(OutputKmers, cc, "_intersect")
  UnionOutput=paste0(OutputKmers, cc, "_Union")
  
  IntersectOpsFileID = paste0(cc, "_ingroups.txt")
  UnionOpsFileID = paste0(cc, "_outgroups.txt")
  
  ShellScriptOutput = paste0(OutputFolder, cc, "_commands.sh")
  
  UniqueOutput = paste0(OutputKmers, cc, "_Unique")
  UniqueOutputText=paste0(OutputKmers, cc, "_MarkerKmers.txt")
  
  lineListShell = c("#!bin/bash","# Performing KMC operations for clade-specific markers","", "source ~/mambaforge/bin/activate ~/mambaforge/envs/pankmer", "")
  lineListShell = append(lineListShell, paste0("kmc_tools complex ",IntersectOpsFileID ))
  lineListShell = append(lineListShell, paste0("kmc_tools complex ",UnionOpsFileID ))
  lineListShell = append(lineListShell, paste0("kmc_tools simple ",IntersectOutput," ",UnionOutput, " kmers_subtract ",UniqueOutput ))
  lineListShell = append(lineListShell, paste0("rm ",IntersectOutput,"*" ))
  lineListShell = append(lineListShell, paste0("rm ",UnionOutput, "*"))
  lineListShell = append(lineListShell, "")
  lineListShell = append(lineListShell, paste0("kmc_tools transform ", UniqueOutput, " dump ",UniqueOutputText))
  
  writeLines(lineListShell, ShellScriptOutput)
  
  

}

# kmc_tools transform $outputUnique dump $outputText
