# Amy Campbell
# July 2023
# Systematic genomic comparisons 

library(stringr)
library(dplyr)


# Comparisons_patient_141_CC1Staphylokinase_1.csv
InputDir = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Distances/"

completenessInfo = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/SummarySequencingMethod.csv")

bashfilefolder="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/IntraPatientTreeScripts/NucDiffBash/"
core_gene_lists="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/IntraPatientTreeScripts/CoreGeneLists/"
  
GenomeFilePathLPC ="/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/"
GenomeFilePathSuffix="_Final.fasta"
nucDiffLPCoutput="/home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/"
phages=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/CDHit/PhagePresenceAbsence.csv")
plasmids=read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Plasmids/Filtered_Plasmid_Presence_Absence.csv")
genes = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv")
outputDFfolder="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/IntraPatientGeneticComparisons/"
plasmids$X.1 = NULL
genes = genes %>% select(Gene,colnames(genes)[grepl("DORN", colnames(genes))])
JustGenomes = genes[,2:ncol(genes)] 
JustGenomes[JustGenomes==""] <- 0
JustGenomes[JustGenomes!=0] <- 1
JustGenomes = JustGenomes %>% mutate_all(function(x) as.numeric(as.character(x)))
genes[2:ncol(genes)] = JustGenomes
genes$NumGenomes = rowSums(genes[2:ncol(genes)])
Variable = genes %>% filter(NumGenomes <220)
Variable$NumGenomes = NULL

filenamelist = list.files(InputDir)
for(inputfilename in filenamelist){
  comparisonDF = read.csv(paste0(InputDir, inputfilename))
  comparisonDF$phenotype = comparisonDF[,3]
  AllCombos = combn(unique(comparisonDF$Cluster), 2)
  patient=str_split(inputfilename, "_")[[1]][3]
  CC = str_split(inputfilename, "_")[[1]][4]
  Phenotype = str_split(inputfilename, "_")[[1]][5]
  comparisonOutputDF = data.frame()
  for(comboIndex in 1:ncol(AllCombos)){
    linelistShell=c("#!bin/bash", "# Making nucdiff scripts for intrapatient comparisons", "mamba ~/mambaforge/bin/activate NucDiffEnv", "")
    
    clusters=AllCombos[,comboIndex]
    shellscriptfile = paste0(bashfilefolder, paste(c(patient, CC, clusters,Phenotype,"NucDiff.sh" ), collapse="_"))
    DF_file = paste0(outputDFfolder, paste(c(patient, CC, clusters,Phenotype,"_GeneComparisons.csv" ), collapse="_"))
    core_gene_file = paste0(core_gene_lists, paste(c(patient, CC, clusters,Phenotype,"_CoreGenes.txt"), collapse="_"))
    comparisonDF_subset = comparisonDF %>% filter(Cluster %in% clusters)
    
    # Cluster with higher mean phenotype value 
    Highcluster = ((comparisonDF_subset %>% group_by(Cluster) %>% summarize(meanval = mean(phenotype)) %>%
                      filter(meanval == max(meanval))))$Cluster
    Lowcluster = setdiff(clusters, Highcluster)
    
    HighclusterGenomes = sapply(( comparisonDF_subset %>% filter(Cluster==Highcluster))$IsolateID, function(x) str_replace(x,"SA", "DORN"))
    LowclusterGenomes = sapply(( comparisonDF_subset %>% filter(Cluster==Lowcluster))$IsolateID, function(x) str_replace(x,"SA", "DORN"))
    
    HighclusterRep = (completenessInfo %>% filter(Genome %in% HighclusterGenomes) %>% arrange(desc(TypeAssembly)))[1,"Genome"]
    HighclusterRep_method = (completenessInfo %>% filter(Genome %in% HighclusterGenomes) %>% arrange(desc(TypeAssembly)))[1,"TypeAssembly"]
    
    # Any phages differentially present?
    ####################################
    phages_low = phages %>% filter(X %in% LowclusterGenomes) %>% select(-X)
    phages_high = phages %>% filter(X %in% HighclusterGenomes) %>% select(-X)
    
    PhageDFHighLow = data.frame(Name=colnames(phages_low),presentlow=colSums(phages_low), presenthigh=colSums(phages_high) )
    Differential = PhageDFHighLow %>% filter((presentlow==0 & presenthigh==nrow(phages_high)) | (presentlow==nrow(phages_low) & presenthigh==0) )  
    
    if(nrow(Differential)>0){
      Differential$PresentOrAbsent_in_High = if_else(Differential$presenthigh > 0, "Present", "Absent")
      NewRows = data.frame(Patient=patient, CCgroup=CC, phenotype=Phenotype,
                           LowCluster=Lowcluster, HighCluster=Highcluster,PresentOrAbsent_in_High=Differential$PresentOrAbsent_in_High,
                           VariantType="Phage", VariantID=Differential$Name,
                           HighClusterReferenceGenome=HighclusterRep, ReferenceGenomeSequenceMethod=HighclusterRep_method)
      
      comparisonOutputDF = rbind(comparisonOutputDF, NewRows)
      
    }
    
    
    # Any plasmids differentially present? 
    plasmids_low = plasmids %>% filter(X %in% LowclusterGenomes) %>% select(-X)
    plasmids_high = plasmids %>% filter(X %in% HighclusterGenomes) %>% select(-X)
    
    PlasmidsDFHighLow = data.frame(Name=colnames(plasmids_low),presentlow=colSums(plasmids_low), presenthigh=colSums(plasmids_high) )
    DifferentialPlasmids = PlasmidsDFHighLow %>% filter((presentlow==0 & presenthigh==nrow(plasmids_high)) | (presentlow==nrow(plasmids_low) & presenthigh==0) )  
    
    if(nrow(DifferentialPlasmids)>0){
      DifferentialPlasmids$PresentOrAbsent_in_High = if_else(DifferentialPlasmids$presenthigh > 0, "Present", "Absent")
      NewRows = data.frame(Patient=patient, CCgroup=CC, phenotype=Phenotype,
                           LowCluster=Lowcluster, HighCluster=Highcluster,PresentOrAbsent_in_High=DifferentialPlasmids$PresentOrAbsent_in_High,
                           VariantType="Plasmid", VariantID=DifferentialPlasmids$Name,
                           HighClusterReferenceGenome=HighclusterRep, ReferenceGenomeSequenceMethod=HighclusterRep_method)
      comparisonOutputDF = rbind(comparisonOutputDF, NewRows)
    }
    
    # Any genes differentially present?
    genes_low = genes %>% select(Gene,LowclusterGenomes)
    genes_high = genes %>% select(Gene,HighclusterGenomes)
    
    genes_low$numgenomes_low = rowSums(genes_low[2:ncol(genes_low)])
    genes_high$numgenomes_high = rowSums(genes_high[2:ncol(genes_high)])
    GenesHighLow = data.frame(Gene=genes_low$Gene,presentlow=genes_low$numgenomes_low , presenthigh=genes_high$numgenomes_high )
    
    DifferentialGenes = GenesHighLow %>% filter((presentlow==0 & presenthigh==length(HighclusterGenomes)) | (presentlow==length(LowclusterGenomes) & presenthigh==0) )  
    CoreGenes = GenesHighLow %>% filter(presenthigh==length(HighclusterGenomes) & (presentlow==length(LowclusterGenomes)))
    write.csv(CoreGenes$Gene, core_gene_file, quote=F,  row.names=F)
    if(nrow(DifferentialGenes)> 0){
      DifferentialGenes$PresentOrAbsent_in_High = if_else(DifferentialGenes$presenthigh > 0, "Present", "Absent")
      NewRows = data.frame(Patient=patient, CCgroup=CC, phenotype=Phenotype,
                           LowCluster=Lowcluster, HighCluster=Highcluster,PresentOrAbsent_in_High=DifferentialGenes$PresentOrAbsent_in_High,
                           VariantType="Gene", VariantID=DifferentialGenes$Gene,
                           HighClusterReferenceGenome=HighclusterRep, ReferenceGenomeSequenceMethod=HighclusterRep_method)
      comparisonOutputDF = rbind(comparisonOutputDF, NewRows)
      
      
    }
    
    # write out bash script for nucdiff
    for(genome in setdiff(c(LowclusterGenomes, HighclusterGenomes), HighclusterRep) ){
      comparestring=paste(HighclusterRep, genome, sep="_")
      commandstring = paste0("nucdiff ", GenomeFilePathLPC, HighclusterRep, GenomeFilePathSuffix,
                             " ",GenomeFilePathLPC, genome, GenomeFilePathSuffix, " ",nucDiffLPCoutput,comparestring, " ", comparestring )
      
      linelistShell= append(linelistShell,commandstring )
      linelistShell= append(linelistShell,"" )
      
    }
    writeLines(linelistShell, shellscriptfile)
    
    write.csv(comparisonOutputDF,file=DF_file )
  }
  
  
}

