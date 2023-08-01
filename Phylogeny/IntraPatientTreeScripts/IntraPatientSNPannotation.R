library(Biostrings)
library(dplyr)
library(stringr)


snpfileprefix = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/AllSNPsites/"
Reffileprefix = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/AllReferences/"
alignmentprefix="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/CoreGeneAlignments/"

snpfilename="extdb_pgaptmp_001954_116_CC5_1_2_Staphylokinase_snpsites"
AnnotatedSNPs = data.frame()
# cols of AnnotatedSNPs:
#     Gene 
#     Phenotype
#     CC
#     Patient
#     Cluster1
#     Cluster2
#     PosStart
#     vartype(INDEL_inFrame, INDEL_FrameShift, AA_Sub, Synonymous, StopCodonSub, INDEL_inFrame_STOP); INDEL_inFrame_STOP refers to an indel where a stop codon is inserted or lost in the indel 
#     PresentOrAbsent_High --> For simple AA substitutions or synonymous mutations, 'AbsentAAsubstitution'/ AbsentSynonymous is the default since I used
#                              a 'high' reference for each comparison
#                              However, for cases where one codon corresponds to "STOP" and the other to an AA, call it a StopCodonSub (stop gained)
#                              and say which (StopGainedHigh if refAA=stop, StopGainedLow if altAA=stop)
#                             Similarly, if there's an indel: if refAllele =* call it 'GapHigh'; if altallele = * call it GapLow
#     Total # consistent SNP/INDEL variants for this comparison (calculate at the end)



snpfilepath = paste0(snpfileprefix, snpfilename)
snpfile = read.table(snpfilepath,skip=3,sep='\t', comment.char="$",header=T)

alignmentfilename = str_replace(snpfilename, "snpsites", "MAFFT.aln")
###########################################
# Parse comparison info based on file name
###########################################
comparisonlist = str_split(snpfilename, "_")[[1]]
genename = paste(comparisonlist[(1: (length(comparisonlist)-6))], collapse="_")
patientID = comparisonlist[length(comparisonlist)-5]
CClabel = comparisonlist[length(comparisonlist)-4]
Cluster1 = comparisonlist[length(comparisonlist)-3]
Cluster2 = comparisonlist[length(comparisonlist)-2]
Phenotype = comparisonlist[length(comparisonlist)-1]


ReferenceFileName = paste0(paste(c("Reference", genename,  patientID,CClabel, Cluster1, Cluster2, Phenotype), collapse="_"), ".fasta")
ReferenceFilePath= paste0(Reffileprefix, ReferenceFileName)
############
# INDELS
############
# First, check for stretches of adjacent insertions/deletions that are identically distributed
deletiongroups =c()
for(j in 1:3500){ # longest length of any snp-sites file is 3097
  deletiongroups[[j]] = c(0)
}
deletionindex=0

# Filter multi-allelic sites since they, by definition, won't be all one allele in low, all another allele in high 
snpfile = snpfile %>% filter(!grepl(",",ALT))

snpfileHigh = snpfile %>% select(colnames(snpfile)[grepl("_High",colnames(snpfile))])
snpfileLow = snpfile %>% select(colnames(snpfile)[grepl("_Low",colnames(snpfile))])

consistent = intersect(which(all(snpfileHigh==0)), which(all(snpfileLow==1)))

consistentSNPs = snpfile[consistent,]

RemainingSNPs = consistentSNPs
for(i in 1:(nrow(consistentSNPs))){
  
  
  if(consistentSNPs[i,"ALT"] == "*"){
    if(i==1){
      deletionindex = deletionindex+1
      deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
    }
    else{
      if((consistentSNPs[i-1, "ALT"] != "*")){
        deletionindex= deletionindex +  1
        deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
      }else{
        # Add to the existing stretch of indels if it's adjacent to its neighbor(which is also a deletion)
        if( ((consistentSNPs[i,"POS"] - consistentSNPs[i-1, "POS"])==1 ) & (consistentSNPs[i-1, "ALT"] == "*") & ( all((consistentSNPs[i,11:ncol(consistentSNPs)]) == (consistentSNPs[i-1,11:ncol(consistentSNPs)])) )  ) {
          deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
        }else{
          deletionindex= deletionindex +  1 
          deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
        }
      }
    }
  }
  
  if(consistentSNPs[i,"REF"] == "*"){
    if( (i==1)){
      deletionindex = deletionindex+1
      deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
      
    }else{
      if(consistentSNPs[i-1, "REF"] != "*"){
        deletionindex= deletionindex +  1
        deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
      }else{
        if( ((consistentSNPs[i,"POS"] - consistentSNPs[i-1, "POS"])==1 ) & (consistentSNPs[i-1, "REF"] == "*") & ( all((consistentSNPs[i,11:ncol(consistentSNPs)]) == (consistentSNPs[i-1,11:ncol(consistentSNPs)])) )  ) {
          deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
        }else{
          deletionindex= deletionindex +  1 
          deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
        }
      }
      
      # Add to the existing stretch of indels if it's adjacent to its neighbor(which is also a deletion)
   
    }
    
  }
}


deletiongroups = lapply(deletiongroups,function(x) setdiff(x, 0) )
deletiongroups = deletiongroups[!isEmpty(deletiongroups)]





for(l in deletiongroups) {
  
    start = min(l)
    stop = max(l)
    lengthIndel = (stop-start) + 1
    
    # Check to see if it's inframe. if it is, check whether there's a 
    # stop codon inserted/deleted
    if( (lengthIndel %% 3) == 0 ){

      
      # check to see if there's a stop codon inserted/deleted
      rows_indel = consistentSNPs %>% filter(POS %in% l)
      rows_indelRef = paste(rows_indel$REF, collapse="")
      rows_indelAlt = paste(rows_indel$ALT, collapse="")
      
      if(grepl("\*", rows_indelAlt)){
        Refcodons = sapply(seq(from=1,to=nchar(rows_indelRef), by=3), function(x) substr(rows_indelRef,x,x+2))
        variableAAs = sapply(Refcodons,function(x) Biostrings::GENETIC_CODE[[x]])
        
      }else{
        Altcodons = sapply(seq(from=1,to=nchar(rows_indelAlt), by=3), function(x) substr(rows_indelAlt,x,x+2))
        variableAAs = sapply(Altcodons,function(x) Biostrings::GENETIC_CODE[[x]])
      }
      if("*" %in% variableAAs){
        vartype="INDEL_inFrame_STOP"
      }else{
        vartype="INDEL_inFrame"
      }
    }else{
      vartype="INDEL_FS"
    }
    
    # If the ALT allele is *, the 'gap' is in the low phenotype;
    # otherwise, it's in the high
    
    testrow = consistentSNPs %>% filter(POS==l[1])
    if(testrow["ALT"]=="*"){
      PresentOrAbsent_High = "GapLow"
    }else{
      PresentOrAbsent_High = "GapHigh"
    }
   
  rowAnnotated = c(genename, Phenotype, CClabel, patientID, Cluster1, Cluster2,start, vartype, PresentOrAbsent_High)  
  AnnotatedSNPs = rbind(AnnotatedSNPs, rowAnnotated)
  
  
  RemainingSNPs = RemainingSNPs %>% filter(!(POS %in% l))
  
  }
  


# Then, if there's anything left in RemainingSNPs after dealing with the indels
###############################################################################

# A SNP position refers to where in the ALIGNMENT a SNP falls; not where in the reference sequence
# itself the snp falls; need to map the reference positions to their alignment positions

if(nrow(RemainingSNPs) > 0){
  
  refidentifier =(names(refSequence))
  
  
  refSequence = Biostrings::readDNAStringSet(ReferenceFilePath, format="fasta")
  alignment = Biostrings::readBStringSet(paste0(alignmentprefix,alignmentfilename), format="fasta")
  alignment = alignment[refidentifier]
  alignmentbases = str_split(as.character(unlist(alignment)), "")[[1]]
  alignmentbases = sapply(alignmentbases, toupper)
  
  indexer =1 
  alignmentpositions = 1:length(alignmentbases)
  refpositions = c()
  for(b in alignmentbases ){
    if(b == "*"){
      refpositions = append(refpositions, "NA")
    }else{
      refpositions = append(refpositions, indexer)
      indexer=indexer+1
      
    }
  }
  
  MapRef_Alignment = data.frame(POS = alignmentpositions,RefPos = refpositions)
  
  
  
  RefBases = str_split(as.character(unlist(refSequence)), "")[[1]]
  AANums =sapply(1:length(RefBases), function(x) ceiling(x/3))
  refDF = data.frame(RefPos = 1:length(RefBases), refbase = RefBases, AAnum =AANums, codonPos = rep(c(1,2,3), length(RefBases)/3))
  refDF$RefPos = sapply(refDF$RefPos, as.character)
  refDF = MapRef_Alignment %>% left_join(refDF,by="RefPos")
  refDF = refDF %>% left_join(RemainingSNPs %>% select(POS, REF, ALT), by="POS")
  row.names(refDF) = refDF$POS
  
  
  for(snpInd in 1:nrow(RemainingSNPs)){
    position = RemainingSNPs[snpInd, "POS"]
    rowinfo = refDF[position,]
    codonbases = (refDF %>% filter(AAnum == as.character(rowinfo["AAnum"])))
    referencecodonbases = paste(codonbases$refbase, collapse="")
    variantcodonbases  = paste(if_else(is.na(codonbases$ALT), codonbases$refbase, codonbases$ALT), collapse="")
    print(referencecodonbases)
    old = Biostrings::GENETIC_CODE[[referencecodonbases]]
    print(variantcodonbases)
    
    new = Biostrings::GENETIC_CODE[[variantcodonbases]]
    presabs = "absent"
    
    if(old==new){
      vartype = "Synonymous"
    }else{
      if(new=="*"){
        vartype = "AA_sub_STOPgained"
        
      }else{
        if(old=="*"){
          vartype = "AA_sub_STOPgained"
          presabs = "present"
        }else{
          vartype = "AA_sub"
          
        }
      }
    }
    rowAnnotated = c(genename, Phenotype, CClabel, patientID, Cluster1, Cluster2,position, vartype, presabs)  
    AnnotatedSNPs = rbind(AnnotatedSNPs, rowAnnotated)
  }
  }
  
  

if(nrow(AnnotatedSNPs) > 0 ){
  AnnotatedSNPs$NumTotalSeqVariants = nrow(AnnotatedSNPs)
}

