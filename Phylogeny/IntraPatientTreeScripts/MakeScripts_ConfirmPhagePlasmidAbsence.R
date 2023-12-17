# Amy Campbell
# 2023
# Identifying which phage/plasmids that were 'variably present/absent' within 
# groups of sequentially cultured isolates from the same lineage/patient
library(stringr)
library(dplyr)
library(seqinr)

phagefastapath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/Hclust_Multifastas/"
plasmidfastapath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/RepMultifastas/"

PhylogeneticClusters = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/IntraPatientPhylogeneticClusters.csv")

outputscriptprefix="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/IntraPatientTreeScripts/ReadMapHGEs/"
outputreferenceprefix = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/ReferencePhagesPlasmids/"
length(unique(PhylogeneticClusters$Cluster))
lengthoutput = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Lengths_HGEs.txt"

Phage_Presence_Absence =  read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/PresenceAbsenceHClustPhages.csv")
Plasmid_Presence_Absence =  read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/Plasmid_Presence_Absence_UTD.csv")

Phage_Presence_Absence$X = sapply(Phage_Presence_Absence$Genome, function(x) str_replace(x, "DORN", "SA"))
Plasmid_Presence_Absence$X = sapply(Plasmid_Presence_Absence$X , function(x) str_replace(x, "DORN", "SA"))

phageDict = c()
plasDict = c()
PresAbsDF= data.frame()

for (clust in unique(PhylogeneticClusters$Cluster)){
  subset = PhylogeneticClusters %>% filter(Cluster==clust)
  plasmids_variable = Plasmid_Presence_Absence %>% filter(X %in% subset$X)
  colsums_plasmid = colSums(plasmids_variable[,3:ncol(plasmids_variable)])
  if(any(colsums_plasmid != 0 & colsums_plasmid!=nrow(plasmids_variable))){
    variableplasmids=colsums_plasmid[colsums_plasmid != 0 & colsums_plasmid!=nrow(plasmids_variable)]

    
    for(plasmiditem in names(variableplasmids)){
      plascols=plasmids_variable[,c("X", plasmiditem)]
      ordered = plascols[rev(order(plascols[,2])),]
      ordered$variable=plasmiditem
      colnames(ordered) = c("X", "PresenceAbsence", "variable")
      PresAbsDF = rbind(PresAbsDF, ordered)
      if(plasmiditem %in% names(plasDict)){
        print("yeah")
        plasDict[[plasmiditem]] = append(plasDict[[plasmiditem]], paste(ordered$X, collapse="_"))
      }else{
        plasDict[[plasmiditem]] =  c(paste(ordered$X, collapse="_"))
      }
    }
    
    
  }
  
  phage_variable = Phage_Presence_Absence %>% filter(X %in% subset$X)
  colsums_phage = colSums(phage_variable[,3:ncol(phage_variable)])
  if(any(colsums_phage != 0 & colsums_phage!=nrow(phage_variable))){
    variablephage=colsums_phage[colsums_phage != 0 & colsums_phage!=nrow(phage_variable)]
    for(phageitem in names(variablephage)){
      phagecols=phage_variable[,c("X", phageitem)]
      ordered = phagecols[rev(order(phagecols[,2])),]
      ordered$variable=phageitem
      colnames(ordered) = c("X", "PresenceAbsence", "variable")
      PresAbsDF = rbind(PresAbsDF, ordered)
      if(phageitem %in% names(phageDict)){
        phageDict[[phageitem]] = append(phageDict[[phageitem]], paste(ordered$X, collapse="_"))
      }else{
        phageDict[[phageitem]] =  c(paste(ordered$X, collapse="_"))
      }
    }
    
  }
}

lengthlist=c("sequence,length")
for(item in names(phageDict)){
  fpath = (paste0(phagefastapath, item, "_all.fasta"))
  fasta= read.fasta(file = fpath)
  
  names(fasta) = sapply( names(fasta),function(x) str_split(x,"_")[[1]][2])
  names(fasta) = sapply( names(fasta), function(x) str_replace(x, "DORN","SA" ))
  for(subitem in phageDict[[item]]){
    # first genome in the list is the 'reference sequence' 
    referenceGenome = str_split(subitem, "_")[[1]][1]
    sequenceout = toupper(fasta[[referenceGenome]])
    lengthlist = append(lengthlist, paste0(paste(item,referenceGenome, sep="_"), ",", nchar(sequenceout)))
    
    # First, write out the reference sequence
    outputfname=paste0(paste(item,referenceGenome, sep="_"), ".fasta")
    write.fasta(sequences = sequenceout, names=paste(item,referenceGenome, sep="_"), file.out = paste0(outputreferenceprefix,outputfname))
    markerspath="/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/Markers/"
    
    scriptoutput=paste0(outputscriptprefix, paste(item,referenceGenome, sep="_"), "_Align.sh")
    linelistScript = c("source ~/mambaforge/bin/activate BowtieEnv23\n",
                       "bowtiepath=/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/BTdbs/\n",
                       "mkdir -p $bowtiepath\n",
                       "export BOWTIE2_INDEXES=$bowtiepath\n",
                       "\n",
                       paste0("outputcounts=",paste(item,referenceGenome, sep="_"),"_Coverage.txt"),
                       "touch $outputcounts",
                       "trimmed1ext=\"trimmedgalore_val_1.fastq\"\n",
                       "samext=\".sam\"\n",
                       "bamext=\".bam\"\n",
                       "blank=\"\"\n",
                       paste0("IndexName=","\"", paste(item,referenceGenome, sep="_"),"\"", "\n"),
                       paste0("markerpath=", markerspath,outputfname ),
                       "bowtie2-build $markerpath $IndexName",
                       "mv *.bt2 $bowtiepath",
                       "\n"
                       )
    allgenomes=str_split(subitem, "_")[[1]]
    
    for(genome in allgenomes){
      genomeNum=str_remove(genome, "SA")
      readspath1=paste0("/project/grice/storage/DFUShortReads2022/trimmedreads/DORN", genomeNum, "trimmedgalore_val_1.fastq")
      readspath2=paste0("/project/grice/storage/DFUShortReads2022/trimmedreads/DORN", genomeNum, "trimmedgalore_val_1.fastq")
      
      linelistScript=append(linelistScript, paste0("readspath1=",readspath1))
      linelistScript=append(linelistScript, paste0("readspath2=",readspath2))
      linelistScript = c(linelistScript, c("basenamefile=$(basename $readspath1)","noext=${basenamefile/$trimmed1ext/$blank}",
                                           "bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext",
                                           "samtools sort $outputpath$noext$samext > $outputpath$noext$bamext", 
                                           "samtools index $outputpath$noext$bamext\n",
                                           "totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)",
                                           "echo $noext\t$totalBasesCovered10x >> $outputcounts",
                                           "\n",
                                           "rm $outputpath$noext$samext",
                                           "rm $outputpath$noext$bamext",
                                           paste0("rm $outputpath$noext$bamext", ".bai")))
      
      
    }
    writeLines(linelistScript, scriptoutput)
    # Then, make a shell script to map reads against the reference sequence 
    
  }
}




for(item in names(plasDict)){
  
  fpath = (paste0(plasmidfastapath, item, "_sequences.fasta"))
  fasta= read.fasta(file = fpath)
  
  names(fasta) = sapply( names(fasta),function(x) str_split(x,"_")[[1]][1])
  names(fasta) = sapply( names(fasta), function(x) str_replace(x, "DORN","SA" ))
  for(subitem in plasDict[[item]]){
    # first genome in the list is the 'reference sequence' 
    referenceGenome = str_split(subitem, "_")[[1]][1]
    sequenceout = toupper(fasta[[referenceGenome]])
    lengthlist = append(lengthlist, paste0(paste(item,referenceGenome, sep="_"), ",", nchar(sequenceout)))
    # First, write out the reference sequence
    outputfname=paste0(paste(item,referenceGenome, sep="_"), ".fasta")
    write.fasta(sequences = sequenceout, names=paste(item,referenceGenome, sep="_"), file.out = paste0(outputreferenceprefix,outputfname))
    markerspath="/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/Markers/"
    
    scriptoutput=paste0(outputscriptprefix, paste(item,referenceGenome, sep="_"), "_Align.sh")
    linelistScript = c("source ~/mambaforge/bin/activate BowtieEnv23\n",
                       "bowtiepath=/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/BTdbs/\n",
                       "mkdir -p $bowtiepath\n",
                       "export BOWTIE2_INDEXES=$bowtiepath\n",
                       "\n",
                       paste0("outputcounts=",paste(item,referenceGenome, sep="_"),"_Coverage.txt"),
                       "touch $soutputcounts",
                       "trimmed1ext=\"trimmedgalore_val_1.fastq\"\n",
                       "samext=\".sam\"\n",
                       "bamext=\".bam\"\n",
                       "blank=\"\"\n",
                       paste0("IndexName=","\"", paste(item,referenceGenome, sep="_"),"\"", "\n"),
                       paste0("markerpath=", markerspath,outputfname ),
                       "bowtie2-build $markerpath $IndexName",
                       "mv *.bt2 $bowtiepath",
                       "\n"
    )
    allgenomes=str_split(subitem, "_")[[1]]
    
    for(genome in allgenomes){
      genomeNum=str_remove(genome, "SA")
      readspath1=paste0("/project/grice/storage/DFUShortReads2022/trimmedreads/DORN", genomeNum, "trimmedgalore_val_1.fastq")
      readspath2=paste0("/project/grice/storage/DFUShortReads2022/trimmedreads/DORN", genomeNum, "trimmedgalore_val_1.fastq")
      
      linelistScript=append(linelistScript, paste0("readspath1=",readspath1))
      linelistScript=append(linelistScript, paste0("readspath2=",readspath2))
      linelistScript = c(linelistScript, c("basenamefile=$(basename $readspath1)","noext=${basenamefile/$trimmed1ext/$blank}",
                                           "bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext",
                                           "samtools sort $outputpath$noext$samext > $outputpath$noext$bamext", 
                                           "samtools index $outputpath$noext$bamext\n",
                                           "totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)",
                                           "echo $noext\t$totalBasesCovered10x >> $outputcounts",
                                           "\n","rm $outputpath$noext$samext",
                                           "rm $outputpath$noext$bamext",
                                           paste0("rm $outputpath$noext$bamext", ".bai")))
      
      
    }
    writeLines(linelistScript, scriptoutput)
    # Then, make a shell script to map reads against the reference sequence 
    
  }
}

writeLines(lengthlist,lengthoutputs )