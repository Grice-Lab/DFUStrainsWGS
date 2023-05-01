# Amy Campbell
# April 2023
# annotated VCF output (SnpEff) for sigB operon
library(dplyr)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)

# takes in: 
##############
# Path to a given gene's annotated VCF from early timepoint
EarlyTimepointPath=args[1]

# Path to a given gene's annotated VCF from late timepoint 
LateTimepointPath=args[2]

# Path to the basecounts file  for the early timepoint (e.g. kalan01_DFUwgs_147_1.txt)
BaseCounts_earlyPath= args[3]

# Path to the basecounts file for the late timepoint
BaseCounts_latePath= args[4]

# Path to output folder
OutputFolder= args[5]

# String containing the transcript ID of interest (e.g., SAUSA300_2024 for rsbV) 
TranscriptString= args[6]

# # Temporary
# ###########
# EarlyTimepointPath = "/Users/amycampbell/Documents/DataInputGithub/data/SigBVariantsMetagenomes/AnnotatedSNPs/kalan01_DFUwgs_106_1_LAC_rsbU.vcf"
# LateTimepointPath = "/Users/amycampbell/Documents/DataInputGithub/data/SigBVariantsMetagenomes/AnnotatedSNPs/kalan01_DFUwgs_106_6_LAC_rsbU.vcf"
# BaseCounts_earlyPath = "/Users/amycampbell/Documents/DataInputGithub/data/SigBVariantsMetagenomes/AnnotatedSNPs/kalan01_DFUwgs_106_1.txt"
# BaseCounts_latePath = "/Users/amycampbell/Documents/DataInputGithub/data/SigBVariantsMetagenomes/AnnotatedSNPs/kalan01_DFUwgs_106_6.txt"
# OutputFolder = "/Users/amycampbell/Documents/DataInputGithub/data/SigBVariantsMetagenomes"
# Example call:
# Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_141_1_LAC_rsbU.vcf
# Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_141_6_LAC_rsbU.vcf
# BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_141_1.txt
# BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_141_6.txt
# Output=/home/acampbe/DFU/data/StressOperonMarkers/VariantComparisons
# Transcriptstring="SAUSA300_2025"
# Rscript SnpEffOutput.R $Early $Late $BaseCounts_earlyPath $BaseCounts_latePath $Output $Transcriptstring
# ############


EarlyTimepoint = read.table(EarlyTimepointPath)
LateTimepoint = read.table(LateTimepointPath)

BaseCounts_early = read.table(BaseCounts_earlyPath) %>% filter(V1=="Chromosome")
BaseCounts_late = read.table(BaseCounts_latePath) %>% filter(V1=="Chromosome")
TranscriptString="SAUSA300_2025"


EarlyTPnum = str_split(basename(EarlyTimepointPath),"_")[[1]][4]
LateTPnum = str_split(basename(LateTimepointPath),"_")[[1]][4]

PatientNum = str_split(basename(EarlyTimepointPath),"_")[[1]][3]
GeneName =  str_remove(str_split(basename(EarlyTimepointPath),"_")[[1]][6], ".vcf")

OutputPrefix = paste0(GeneName, "_Variants_",PatientNum, "_", EarlyTPnum, "_",LateTPnum )
OutputFiltered = paste0(OutputPrefix, "_Final.tsv")
OutputLowDepth = paste0(OutputPrefix, "_LowDepth.tsv")

OutputPathFiltered = file.path(OutputFolder, OutputFiltered)
OutputPathLowDepth = file.path(OutputFolder, OutputLowDepth)


EarlyTimepoint$VariantIdentifier = paste(EarlyTimepoint$V2, EarlyTimepoint$V5, sep="_")
LateTimepoint$VariantIdentifier = paste(LateTimepoint$V2, LateTimepoint$V5, sep="_")

SharedVars = intersect(EarlyTimepoint$VariantIdentifier, LateTimepoint$VariantIdentifier)

# Filter out variants present in both
EarlyTimepoint = EarlyTimepoint %>% filter(!(VariantIdentifier %in% SharedVars))
LateTimepoint = LateTimepoint %>% filter(!(VariantIdentifier %in% SharedVars))


#####################################################################################
# Parse the two timepoints' VCF files
#####################################################################################
EarlyTimepoint$DP = sapply(EarlyTimepoint$V8, function(x) str_split(x, ";")[[1]][1])
EarlyTimepoint$I16 = sapply(EarlyTimepoint$V8, function(x) str_split(x, ";")[[1]][2])
EarlyTimepoint$ANN = sapply(EarlyTimepoint$V8, function(x) (str_split(x, ";")[[1]])[str_detect((str_split(x, ";")[[1]]), "ANN")]  )
EarlyTimepoint$VariantType = sapply(EarlyTimepoint$ANN, function(x) str_split(x, "\\|")[[1]][2])
EarlyTimepoint$PredictedEffect = sapply(EarlyTimepoint$ANN, function(x) str_split(x, "\\|")[[1]][3])
EarlyTimepoint$Gene = sapply(EarlyTimepoint$ANN, function(x) str_split(x, "\\|")[[1]][5])
EarlyTimepoint$TimepointPresent = EarlyTPnum
EarlyTimepoint$TimepointAbsent = LateTPnum


EarlyTimepoint = EarlyTimepoint %>% select(V2, V4, V5,DP, I16, VariantType, PredictedEffect, Gene, TimepointPresent,TimepointAbsent) %>% filter(Gene==TranscriptString)
colnames(EarlyTimepoint) = c("BasePosition", "Ref", "Alt", "DP", "I16", "VariantType", "PredictedEffect", "Gene", "TimepointPresent", "TimepointAbsent")
  
LateTimepoint$DP = sapply(LateTimepoint$V8, function(x) str_split(x, ";")[[1]][1])
LateTimepoint$I16 = sapply(LateTimepoint$V8, function(x) str_split(x, ";")[[1]][2])
LateTimepoint$ANN = sapply(LateTimepoint$V8, function(x) (str_split(x, ";")[[1]])[str_detect((str_split(x, ";")[[1]]), "ANN")]  )
LateTimepoint$VariantType = sapply(LateTimepoint$ANN, function(x) str_split(x, "\\|")[[1]][2])
LateTimepoint$PredictedEffect = sapply(LateTimepoint$ANN, function(x) str_split(x, "\\|")[[1]][3])
LateTimepoint$Gene = sapply(LateTimepoint$ANN, function(x) str_split(x, "\\|")[[1]][5])
LateTimepoint$TimepointPresent = LateTPnum
LateTimepoint$TimepointAbsent = EarlyTPnum

LateTimepoint = LateTimepoint %>% select(V2, V4, V5,DP, I16, VariantType, PredictedEffect, Gene, TimepointPresent, TimepointAbsent) %>% filter(Gene==TranscriptString)
colnames(LateTimepoint) = c("BasePosition", "Ref", "Alt", "DP", "I16", "VariantType", "PredictedEffect", "Gene", "TimepointPresent", "TimepointAbsent")

FinalDF = rbind(EarlyTimepoint, LateTimepoint)


DepthsEarly = BaseCounts_early %>% filter(V3>=5)
DepthsLate = BaseCounts_late %>% filter(V3>=5)
IncludedBases = intersect(DepthsEarly$V2, DepthsLate$V2 )


OutputFilteredDF = FinalDF %>% filter(BasePosition %in% IncludedBases)
OutputLowDepthDF = FinalDF %>% filter(!(BasePosition %in% IncludedBases))


if(nrow(OutputFilteredDF>0)){
  write.table(OutputFilteredDF, file=OutputPathFiltered)
}else{
  print(paste0("No variants in the ", OutputPrefix, " comparison with >=5X depth in both samples."))
}

if(nrow(OutputFilteredDF>0)){
  write.table(OutputFilteredDF, file=OutputLowDepthDF)
}else{
  print(paste0("No variants in the ", OutputPrefix, " comparison with <5X depth in both samples."))
}
