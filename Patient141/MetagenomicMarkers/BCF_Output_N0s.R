# Amy Campbell 
# Take the bcf output of the SNP calls of reads to markers to build a big .csv

library(dplyr)
library(stringr)


MarkerVals = function(BCF, MarkerPosition, MarkerName){
  
  BCF_DF = BCF %>% filter(POS==MarkerPosition & CHROM==MarkerName)
  if(nrow(BCF_DF) == 0){
    
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

outputfilepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/ParameterRange_N0s.csv"

FDirectory="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BTalignments/TestRangeParameters/N0Tests/"
AllFiles = list.files(FDirectory) 
BCF_files= AllFiles[grepl(AllFiles, pattern=".bcf")]


sourcelist=c()
Markerlist=c()
TotalDepthlist=c()
RefDepthlist=c()
AltDepthlist=c()
dvaluelist=c()
rvaluelist=c()
lvaluelist=c()
qvaluelist=c()
s2list=c()
i=0
for(FilePath in BCF_files){

  filename=basename(FilePath)
  
  # StaphEpireads_30_5_20_75_N0_13.bcf
  
  BCFPrefix <- (str_split(filename, pattern="reads_"))[[1]][1]
  BCFSuffix <- str_remove((str_split(filename, pattern="reads_"))[[1]][2], pattern=".bcf")
  listparams=str_split(BCFSuffix, pattern="_")
  dvalue <- as.numeric(listparams[[1]][1])
  rvalue <-  as.numeric(listparams[[1]][2])
  lvalue <-  as.numeric(listparams[[1]][3])
  qvalue <- as.numeric(listparams[[1]][6])
  s2value <-as.numeric(listparams[[1]][4])
  
  for(marker in c("Marker1", "Marker2", "Marker3")){
    i=i+1
    
    dvaluelist[i] = dvalue
    rvaluelist[i] = rvalue
    lvaluelist[i] = lvalue
    qvaluelist[i] = qvalue
    Markerlist[i] = marker
    sourcelist[i] = BCFPrefix
    s2list[i] = s2value
    
  
  # Get Full path for this file
  FullPath=paste0(FDirectory, FilePath)
  
  # Try opening this file
  trystatement = try(read.delim(FullPath, comment.char = "#", header=F), silent=T)

  if(inherits(trystatement, "try-error")){
    
    TotalDepthlist[i] = 0
    RefDepthlist[i] = 0
    AltDepthlist[i] = 0 
  
    print("Empty File, sRy!")} else{
      
      FileObject=read.delim(FullPath,comment.char = "#", header=F)
      colnames(FileObject) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT", "V10")

      markerlist = MarkerVals(FileObject, 250, marker)
      
      TotalDepthlist[i] = markerlist[2]
      RefDepthlist[i] = markerlist[3]
      AltDepthlist[i] = markerlist[4]
      
      
    }
  }
  
}

OutputDF = data.frame(Source=sourcelist,
                      MarkerID=Markerlist, 
                      TotalUseableDepth=TotalDepthlist,
                      RefAlleleDepth=RefDepthlist,
                      AltAlleleDepth=AltDepthlist,
                      Dvalue=dvaluelist,
                      Rvalue=rvaluelist,
                      Lvalue=lvaluelist,
                      Qvalue=qvaluelist,
                      S2value=s2list)
OutputDF$Nvalue = "N0"
write.csv(OutputDF, file=outputfilepath, quote = F, row.names = F)











