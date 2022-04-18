# Amy Campbell 
# Take the bcf output of the SNP calls of reads to markers to build a big .csv

library(dplyr)
library(stringr)

# function which makes a list of MarkerID, TotalUseableDepth, RefAlleleDepth, AltAlleleDepth

MarkerVals = function(BCF, MarkerPosition, MarkerName){
  
  BCF_DF = BCF %>% filter(POS==MarkerPosition & CHROM==MarkerName)
  if(nrow(BCF_DF) == 0){
    
    # c(MarkerID, TotalUseableDepth, RefAlleleDepth, AltAlleleDepth)
    return(c(MarkerName, 0,0,0))
  }else{
    BCF_DF$InfoDP = sapply(BCF_DF$INFO, function(x) str_split(string=x, pattern = ";")[[1]][2])
    BCF_DF$InfoDP = sapply(BCF_DF$InfoDP, function(x) str_split(string=x, pattern = "I16=")[[1]][2])
    BCF_DF$DepthRef = sapply(BCF_DF$InfoDP, function(x) sum(as.numeric(as.character(str_split(string=x, pattern = ",")[[1]][1:2]))))
    BCF_DF$DepthAlt = sapply(BCF_DF$InfoDP, function(x) sum(as.numeric(as.character(str_split(string=x, pattern = ",")[[1]][3:4]))))
    BCF_DF$InfoDP = BCF_DF$DepthRef+BCF_DF$DepthAlt
    return(c( MarkerName, BCF_DF$InfoDP,BCF_DF$DepthRef, BCF_DF$DepthAlt))  
  }
  
  
}

outputfilepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/ParameterRange.csv"

FDirectory="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BTalignments/TestRangeParameters/"
AllFiles = list.files(FDirectory) 
BCF_files= AllFiles[grepl(AllFiles, pattern=".bcf")]



OutputDF = data.frame()

for(FilePath in BCF_files){
  
  # Get info for this file
  BCFPrefix <- (str_split(filename, pattern="reads_"))[[1]][1]
  BCFSuffix <- str_remove((str_split(filename, pattern="reads_"))[[1]][2], pattern=".bcf")
  listparams=str_split(BCFSuffix, pattern="_")
  dvalue <- as.numeric(listparams[[1]][1])
  rvalue <-  as.numeric(listparams[[1]][2])
  lvalue <-  as.numeric(listparams[[1]][3])
  qvalue <-  as.numeric(listparams[[1]][4])
  
  # Get Full path for this file
  FullPath=paste0(FDirectory, FilePath)
  
  # Try opening this file
  trystatement = try(read.delim(FullPath, comment.char = "#", header=F), silent=T)
  
  # If it's empty other than comment lines, 
  #colnames(OutputDF) = c("Source", "MarkerID", "TotalUseableDepth", "RefAlleleDepth", "AltAlleleDepth", "dvalue", "rvalue", "lvalue", "qvalue")
  
  if(inherits(trystatement, "try-error")){
    
    marker1 = c("Marker1", 0,0,0,dvalue,rvalue,lvalue,qvalue)
    NewRow = append(BCFPrefix,marker1)
    NewRow = append(NewRow, c(dvalue,rvalue,lvalue,qvalue))
    OutputDF = rbind(OutputDF, NewRow)
    
    marker2= c("Marker2", 0,0,0,dvalue,rvalue,lvalue,qvalue)
    NewRow = append(BCFPrefix,marker2)
    OutputDF = rbind(OutputDF, NewRow)
    
    marker3=c("Marker3", 0,0,0,dvalue,rvalue,lvalue,qvalue)
    NewRow = append(BCFPrefix,marker3)
    OutputDF = rbind(OutputDF, NewRow)
    
    print("Empty File, sRy!")} else{
      
      FileObject=read.delim(FullPath,comment.char = "#", header=F)
      colnames(FileObject) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT", "V10")
      filename=basename(FilePath)
      
      marker1 =MarkerVals(FileObject, 250, "Marker1")
      NewRow = append(BCFPrefix,marker1)
      OutputDF = rbind(OutputDF, NewRow)
      
      marker2 =MarkerVals(FileObject, 250, "Marker2")
      NewRow = append(BCFPrefix,marker2)
      OutputDF = rbind(OutputDF, NewRow)
      
      marker3 =MarkerVals(FileObject, 250, "Marker3")
      NewRow = append(BCFPrefix,marker3)
      OutputDF = rbind(OutputDF, NewRow)
      
      
    }
  
}

colnames(OutputDF) = c("Source", "MarkerID", "TotalUseableDepth", "RefAlleleDepth", "AltAlleleDepth", "dvalue", "rvalue", "lvalue", "qvalue")
write.csv(OutputDF, file=outputfilepath, quote = F, row.names = F)











