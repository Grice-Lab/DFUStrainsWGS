# Amy Campbell
# March 2023

# Takes in:
# 1) Path to folder containing the VCFs 
# e.g. 

# 2) case control csv which was input into whatsGNU
# e.g. /Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscript/case_controls_patient176.csv

# 3) Name of the reference genome
# e.g. DORN1863

library(stringr)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)

# takes in: 
##############
# ReferenceGenome="DORN1863"
ReferenceGenome=args[1]

# CaseControlPath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/case_controls_patient176.csv"
CaseControlPath=args[2]

# DirectoryPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/GeneFastas/MAFFT"
DirectoryPath=args[3]

# OutputPath = /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/Patient176_Variants.tsv
OutputPath= args[4]

ReadCaseCtrls=read.csv(CaseControlPath, header=F)
ReadCaseCtrls$clean = sapply(ReadCaseCtrls$V1, function(x) (str_split(x,"_WhatsGNU_report.txt"))[[1]][1] )

Caselist = (ReadCaseCtrls %>% filter(V2=="case"))$clean
Controllist = (ReadCaseCtrls %>% filter(V2=="control"))$clean

FilesInput = list.files(DirectoryPath)
FilesInput = FilesInput[grepl("SNPs_",FilesInput )]


FullVariantDF=data.frame()
FinalGeneList =c()
for(f in FilesInput){
  
  # path to VCF
  fullpath=file.path(DirectoryPath, f)
  
  # name of gene
  InputOrthologName=str_replace(f, "SNPs_", "")
  
  # Read in the VCF
  inputDF = read.table(fullpath,skip=3,sep='\t', comment.char="$",header=T)
  
  # Fixing stupid read table assumption that string T=true lol 
  inputDF$REF = if_else(inputDF$REF==TRUE, "T", inputDF$REF)
  inputDF$ALT = if_else(inputDF$ALT==TRUE, "T", inputDF$ALT)
  
  # Remove the extraneous instance of the reference gene which we used for aligning all of them 
  inputDF = inputDF %>% select(-paste0(ReferenceGenome, ".1"))
  
  # Get the length of the alignment
  inputstring=file(fullpath, open="r")
  inputstringLength=(readLines(inputstring))[2]
  lengthAlign = (str_split(inputstringLength,"length=" ))[[1]][2]
  lengthAlign = as.integer(str_split(lengthAlign, ">")[[1]][1])
  
  # If more than 10% of the alignment's length is variants, I anticipate potential mapping issues 

  if(nrow(inputDF) > .1*lengthAlign){
    print("Too many mismatches. Remove this gene :(")
  } else{
    inputDF = inputDF %>% filter(nchar(ALT)==1)
    if(nrow(inputDF>0)){
      FinalGeneList = append(FinalGeneList,InputOrthologName )
    }
    # Indices at which all 'cases' == 0
    CaseCols = inputDF %>% select(Caselist)
    Req1 = which(rowSums(CaseCols)==0)
    
    # Indices at which all 'controls' == 1
    CtrlCols = inputDF %>% select(Controllist)
    Req2 = which(rowSums(CtrlCols)==ncol(CtrlCols))
    
    # Rows at which all controls == 1 and all cases == 0 
    RowsUse = inputDF[intersect(Req1, Req2),]
    
    # Now we have a little dataframe of the good rows
    DBAdd = RowsUse %>% select(REF, ALT, POS)
    print(DBAdd)
    DBAdd$Gene = InputOrthologName
    colnames(DBAdd) = c("Case", "Control", "Position", "Gene")

    FullVariantDF = rbind(FullVariantDF, DBAdd)
  }
  
}

OutputGeneList = str_replace(OutputPath, ".tsv", "_FinalGeneList.txt")

write.table(FullVariantDF, file=OutputPath, quote=F)

write.table(FullVariantDF, file=OutputPath, quote=F,col.names = F, row.names=F)





      