# Amy Campbell
# 3/26/2023
# Output of bowtie2/samtools/bcftools alignment and variant calling of a 
library(dplyr)
library(stringr)

#BCFPath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/ShotgunScripts/BCFs"
#VariantList="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/ShotgunScripts/Patient176_Variants.tsv"
#PatientID="176"
#CaseString="CC5"
#ControlString="CC1"
# 

# Example run : 
################
# BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient176/output/BCFs"
# VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/Patient176_Variants.tsv
# Patient="176"
# Case="CC5"
# Control="CC1"
# OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions
#
# Rscript $BCFs $VariantListPath $Patient $Case $Control $OutputFolder

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 6){
  print("Not enough arguments(need 6) >:( ")
  exit()
}else{
  BCFPath=args[1]
  
  # Absolute path to the TSV output by ReadVCFs.R listing the variants that are present in controls, absent from cases, when you use a case reference
  VariantList=args[2]
  
  # Numeric ID of patient (e.g., "176")
  PatientID=args[3]
  
  # What clade does 'case' represent
  CaseString=args[4]
  
  # What clade does 'control' represent
  ControlString=args[5]
  
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


VariantsToCheck = read.table(VariantList,header=T)
VariantsToCheck$GenePos= paste(VariantsToCheck$Gene,VariantsToCheck$Position, sep="_" )
VariantsToCheck$GenePosAlt = paste(VariantsToCheck$GenePos, VariantsToCheck$Control, sep="_")


ListOutputs = list.files(BCFPath)

OutputDF=data.frame()
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
  BCF_markers_Final$NumCase = sapply(BCF_markers_Final$DP4Col, function(x) AlleleDepth(x, "Ref") )
  BCF_markers_Final$NumCtrl = sapply(BCF_markers_Final$DP4Col, function(x) AlleleDepth(x, "Alt") )
  TotalMarkersUsed=nrow(BCF_markers_Final)
  TotalControl=sum(BCF_markers_Final$NumCtrl)
  TotalCase=sum(BCF_markers_Final$NumCase)
  if( (TotalControl+TotalCase) == 0){
    print("No Mapping to marker varaints. Skipping this timepoint")
  }else{
    OutputDF= rbind(OutputDF, c(Timepoint, TotalControl,TotalCase, (100*TotalCase/(TotalCase+TotalControl)) ))
    
  }
}

colnames(OutputDF) = c("Timepoint", "ControlBases", "CaseBases", "PctCase")

OutputDF$CaseClade=CaseString
OutputDF$ControlClade=ControlString

OutputFname=paste0("Patient_", PatientID, "_composition.csv")
OutputFile=file.path(OutputFolder,OutputFname)
write.csv(OutputDF, file=OutputFile)

