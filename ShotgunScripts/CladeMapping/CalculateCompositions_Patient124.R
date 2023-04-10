# Amy Campbell
# 3/26/2023
# Output of bowtie2/samtools/bcftools alignment and variant calling of a 
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4){
  print("Not enough arguments(need 4) >:( ")
  exit()
}else{ 
  
  BCFPath=args[1]
  
  # Absolute path to the TSV output by ReadVCFs.R listing the variants that are present in controls, absent from cases, when you use a case reference
  VariantList=args[2]
  
  # Numeric ID of patient (e.g., "176")
  PatientID=args[3]

  
  OutputFolder=args[4]  
}

AlleleDepth = function(DP4String, RefOrAlt){
  # DP4 string is the DP4Col value
  #RefOrAlt can be "Ref" or "Alt" (lol)
  # if Ref, returns total # high quality reference alleles  
  # if Alt, returns total # high quality alterate alleles
  JustDepths = (str_split(DP4String, ":"))[[1]][2]
  DepthsList= str_split(JustDepths, ",")
  DepthsList = sapply(DepthsList, function(x) as.numeric(as.character(x)))
  if(RefOrAlt=="Alt"){
    return(DepthsList[3]+DepthsList[4] )
  }else if(RefOrAlt=="Ref"){
    return(DepthsList[1] + DepthsList[2])
  }else{
    print("Error:  preference of Ref or Alt given improperly")
  }
  
}

ListOutputs = list.files(BCFPath)

OutputDF=data.frame()
VariantsToCheck = read.table(VariantList,header=F)
colnames(VariantsToCheck) = c("RefAllele", "AltAllele", "Position", "Gene","Type")

VariantsToCheck$GenePos= paste(VariantsToCheck$Gene,VariantsToCheck$Position, sep="_" )
VariantsToCheck$GenePosAlt = paste(VariantsToCheck$GenePos, VariantsToCheck$AltAllele, sep="_")


for(bcf in ListOutputs){


  
  
  FirstString=paste0("kalan01_DFUwgs_", PatientID, "_")
  Timepoint=(str_split(string = bcf, pattern = FirstString))[[1]][2]
  Timepoint=(str_split(Timepoint, pattern="_marker.bcf"))[[1]][1]
  
  # Read in the bcf file
  Fpath=file.path(BCFPath, bcf)
  BCFoutput=read.delim(Fpath,comment.char = "#", header=F)
  colnames(BCFoutput) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT", "DP4Col")
  
  # make a variable that's the gene_position format
  BCFoutput$GenePos=paste(BCFoutput$CHROM, BCFoutput$POS, sep="_")
  BCF_markers = BCFoutput %>% filter(GenePos %in% VariantsToCheck$GenePos)
  BCF_markers$alternate = sapply(BCF_markers$ALT, function(x) as.factor(str_split(x, ",")[[1]][1]))
  BCF_markers$GeneposAlt = paste(BCF_markers$GenePos, BCF_markers$alternate, sep="_")
  
  BCF_markers_Final = BCF_markers %>% filter(((GeneposAlt %in% VariantsToCheck$GenePosAlt) | alternate=="<*>" ))
  BCF_markers_Final$NumRef = sapply(BCF_markers_Final$DP4Col, function(x) AlleleDepth(x, "Ref") )
  BCF_markers_Final$NumAlt = sapply(BCF_markers_Final$DP4Col, function(x) AlleleDepth(x, "Alt") )
  TotalMarkersUsed=nrow(BCF_markers_Final)
  
  
  BCF_markers_Final_Typed = BCF_markers_Final 
  BCF_markers_Final_Typed$GenePosAlt = BCF_markers_Final_Typed$GeneposAlt
  BCF_markers_Final_Typed = BCF_markers_Final_Typed %>% left_join(VariantsToCheck, by="GenePos")
  TotalRefAlleles = sum(BCF_markers_Final_Typed$NumRef)
  
  A_cases = BCF_markers_Final_Typed %>% filter(Type=="A")
  B_cases = BCF_markers_Final_Typed %>% filter(Type=="B")
  
  NormalizedCC5 = TotalRefAlleles/(sum(BCF_markers_Final_Typed$NumRef) + sum(BCF_markers_Final_Typed$NumAlt))
  NormalizedCC9 = (sum(A_cases$NumAlt))/(sum(A_cases$NumAlt) + sum(A_cases$NumRef))
  NormalizedCC8 = (sum(B_cases$NumAlt))/(sum(B_cases$NumAlt) + sum(B_cases$NumRef))
  
  PropCC5=NormalizedCC5/(NormalizedCC5+NormalizedCC9+NormalizedCC8)
  PropCC9=NormalizedCC9/(NormalizedCC5+NormalizedCC9+NormalizedCC8)
  PropCC8=NormalizedCC8/(NormalizedCC5+NormalizedCC9+NormalizedCC8)
  
  OutputDF = rbind(OutputDF, c(Timepoint,100*PropCC5, 100*PropCC9,100*PropCC8, (sum(BCF_markers_Final_Typed$NumAlt) + sum(BCF_markers_Final_Typed$NumRef))))
  
  
  
}

colnames(OutputDF) = c("Timepoint", "CC5", "CC9", "CC8", "TotalDepth")


OutputFname=paste0("Patient_", PatientID, "_composition.csv")
OutputFile=file.path(OutputFolder,OutputFname)

write.csv(OutputDF, file=OutputFile)

