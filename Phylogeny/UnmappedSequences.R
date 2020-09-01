# Amy Campbell
# This script goes through the output of:
# (1) MiniMap2 pairwise alignments between each isolate and its
# corresponding "closest reference" on the core genome SNP tree (the mappings for which can be found in ClosestReferences.txt)
# The alignments are in PAF format which is just a tsv. Reports:
#     - Total length of sequence in the isolate that didn't align anywhere in the close reference
#     - Maximum length of stretches in the isolate which didn't map anywhere in the close reference 
# (2) plasmidcounts.txt output by Count_Plasmids.sh
#
# (3) /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/Quast_CleanedDORN/transposed_report.tsv SKIP THE FIRST LINE OF THIS ONE 
# 
library(dplyr)

# Filepaths
###########
Alignmentfilepath="/home/acampbe/DFU/data/WGS_2020/MinimapPairwiseAlignments/"
ClosestRefPath="/home/acampbe/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/ClosestReferences.txt"
PlasmidsPath= "/home/acampbe/Club_Grice/EAGenomeAssembly/CoverageAnalyses/plasmidcounts.txt"
AssemblyStatsPath="/project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/Quast_CleanedDORN/transposed_report.tsv"

ClosestRefFrame = read.csv(ClosestRefPath,header=F, sep='\t')
colnames(ClosestRefFrame) = c("DORN", "Reference")

PlasmidsFrame = read.csv(PlasmidsPath, header=F)
colnames(PlasmidsFrame) = c("DORN", "NPlasmids")
#print(PlasmidsFrame)

AssemblyStatsFrame = read.csv(AssemblyStatsPath, sep='\t', header=F, skip = 1)
colnames(AssemblyStatsFrame) = c("AssemblyName", "NumContigs0", "NumContigs1000",
                                 "NumContigs5000", "NumContigs10000", "NumContigs25000",
                                 "NumContigs50000", "Length0", "Length1000", "Length5000",
                                 "Length10000", "Length25000","Length50000", "Ncontigs", "LargestContig",
                                 "TotalAssemblyLength", "GC", "N50", "N75", "L50", "L75", "NsPer100kb")    
AssemblyStatsFrame$DORN = sapply(AssemblyStatsFrame$AssemblyName, function(x) stringr::str_split(x, "_")[[1]][1])
AssemblyStatsFrame  = AssemblyStatsFrame[c("DORN", "Ncontigs", "LargestContig", "TotalAssemblyLength", "GC", "N50","L50")]
print(colnames(AssemblyStatsFrame))

BigFrame = ClosestRefFrame
BigFrame$MaxLength = rep(NA, (dim(BigFrame))[1]) 

BigFrame$TotalLength = rep(NA, (dim(BigFrame))[1])


  
filelist= list.files(Alignmentfilepath, full.names=T)

for(file in filelist){
  fname = basename(file)
  dornname = stringr::str_replace(fname, ".paf", "" )
  AlignmentFile = read.csv(file, sep='\t', header=F)
  colnames(AlignmentFile) = c("QueryName", "QueryLength","QueryStart","QueryEnd", "Sense", "TargetName", "TargetLength", "TargetStart", "TargetEnd", "NumMatches", "MappingLength", "MappingQual", "V13", "V14", "V15", "V16", "V17", "V18")
  Unmapped = AlignmentFile %>% filter(TargetName=="*")
  if(dim(Unmapped)[1] > 0){

 	TotalLength= sum(Unmapped$QueryLength)
  	MaxLength=max(Unmapped$QueryLength)
	} else{
	TotalLength=0
	MaxLength=0
}
  BigFrame[which(BigFrame$DORN == dornname),  "MaxLength"]<- MaxLength
  BigFrame[which(BigFrame$DORN == dornname), "TotalLength"] <- TotalLength
  
}
BigFrame = BigFrame %>% left_join(AssemblyStatsFrame, by="DORN")
BigFrame = BigFrame %>% left_join(PlasmidsFrame, by="DORN")

write.csv(BigFrame, "BigFrame.csv")
