# Amy Campbell 
# 4-2022
# Extract homologous refseqs from a PGAP .gff 'attributes' column 
# from a .gff output by PGAP & reformatted with with the sequences cut off 


# My defaults for 925 
#gffpath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/DORN925_ref.gff"
#outputpathMap="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RefSeqMapping925.txt"
#outputpathRefs="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RefSeqIDs925.txt"


# Read in Arguments from user 
args <- commandArgs()
gffpath=args[1]
outputpathMap=args[2]
outputpathRefs=args[3]

# Function that takes the attributes column and splits it to pieces 
GetRefSeq = function(string, option){
  # option = 1 -- "ID=" field 
  # option = 2 -- "RefSeq field"
  
  
  stringlist = (str_split(string, pattern=";"))[[1]]
  stringid=str_remove(stringlist[1], pattern="ID=")
  
  
  if(option==1){
    stringitem=stringlist[grepl(stringlist, pattern="AA sequence:RefSeq:")]
    if(length(stringitem) >0){
      returnstring = (str_split(stringitem, pattern="RefSeq:"))[[1]][2]
      return(returnstring)
    }else{
      return(paste0("NoRefSeq_", stringid))
    }
  }else(
    return(stringid)
  )
  
}

# Read the GFF file
gff_file= read.delim(gffpath, sep='\t', skip=3, header=F)

attributescol = gff_file$V9

# Get the locus IDs (ID=<locus ID>)
gff_file$V10 = sapply(attributescol, function(x) GetRefSeq(x, 1))

# get  the RefSeq IDs
gff_file$V11 = sapply(attributescol, function(x) GetRefSeq(x, 2))

outputdf = gff_file %>% select(V11, V10)

# A mapping of locus IDs --> RefSeq IDs
write.table(outputdf$V10, file=outputpathRefs, quote=F, row.names=F,col.names = F)

# Just the locus IDs for upload to Uniprot's mapper
write.table(outputdf, file=outputpathMap, quote=F, row.names=F,col.names = F)


