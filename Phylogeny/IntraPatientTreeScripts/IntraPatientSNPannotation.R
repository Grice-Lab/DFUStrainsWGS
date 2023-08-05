library(Biostrings)
library(dplyr)
library(stringr)


snpfileprefix = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/AllSNPsites/"
Reffileprefix = "/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/AllReferences/"
alignmentprefix="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/CoreGeneAlignments/"
OutputSNPannotations="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/SNPannotations/"

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

BigDF = data.frame()
Big_Indel_DF = data.frame()
list.files(snpfileprefix)
for(snpfilename in list.files(snpfileprefix)){
  print(snpfilename)
  AnnotatedSNPs = data.frame()
  
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
  
  snpfileHighSums = rowSums(snpfileHigh)
  snpfileLowSums = rowSums(snpfileLow)
  
  
  
  consistentHigh  = intersect(which(snpfileHighSums==0), which(snpfileLowSums==ncol(snpfileLow)))
  consistentLow = intersect(which(snpfileLowSums==0), which(snpfileHighSums==ncol(snpfileHigh)))
  consistent = unique(consistentHigh, consistentLow)
  consistentSNPs = snpfile[consistent,]
  
  RemainingSNPs = consistentSNPs
  
  if(nrow(consistentSNPs) > 0 ){
    # look for and index deletions
    for(i in 1:(nrow(consistentSNPs))){
      
      
      if(consistentSNPs[i,"ALT"] == "*"){
        
        # If it's the very first one in the df, make a new entry for it
        if(i==1){
          deletionindex = deletionindex+1
          deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
        }else{
          #If it's the first in a possible string of deletions, make a new entry for it
          if((consistentSNPs[i-1, "ALT"] != "*")){
            deletionindex= deletionindex +  1
            deletiongroups[[deletionindex]] = append(deletiongroups[[deletionindex]], consistentSNPs[i, "POS"] )
          }else{
            # Add to the existing stretch of indels if it's adjacent to its neighbor(which is also a deletion)
            if( ((consistentSNPs[i,"POS"] - consistentSNPs[i-1, "POS"])==1 ) & (consistentSNPs[i-1, "ALT"] == "*") & (all((consistentSNPs[i,11:ncol(consistentSNPs)]) == (consistentSNPs[i-1,11:ncol(consistentSNPs)])))){
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
          
        }
        
      }
    }
  
  }

  # Make a map of reference gene position to alignment position
  ##############################################################3
  refSequence = Biostrings::readDNAStringSet(ReferenceFilePath, format="fasta")
  refidentifier =(names(refSequence))
  
  alignment = Biostrings::readBStringSet(paste0(alignmentprefix,alignmentfilename), format="fasta")
  alignment = alignment[refidentifier]
  alignmentbases = str_split(as.character(unlist(alignment)), "")[[1]]
  alignmentbases = sapply(alignmentbases, toupper)
  
  indexer =1 
  alignmentpositions = 1:length(alignmentbases)
  refpositions = c()
  
  #iterate through the bases of the alignment
  for(b in alignmentbases ){
    
    # If the alignment version of the ref genome="*", then there's a gap in it there. 
    # Append "NA" to the ref positions since the base in the reference genome at that point is NA
    # this does deal with things that are missing from the reference genome (and therefore gaps)
    # but what about when something is inserted into the reference genome and missing from the others? 
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
  MapRef_Alignment$RefPos = sapply(MapRef_Alignment$RefPos, as.character)
  
  refDF = MapRef_Alignment %>% left_join(refDF,by="RefPos")
  row.names(refDF) = refDF$POS
  deletiongroups = lapply(deletiongroups,function(x) setdiff(x, 0) )
  deletiongroups = deletiongroups[!isEmpty(deletiongroups)]
  
  # actually file away/annotate the indel, having ID'd the
  # stretches in which each one exists
  
  for(l in deletiongroups){
      start = min(l)
      stop = max(l)
      lengthIndel = (stop-start) + 1
      
      rows_indel_all = consistentSNPs %>% filter(POS %in% l)
      RefBases = paste(rows_indel_all$REF, collapse="")
      AltBases = paste(rows_indel_all$ALT, collapse="")
      position_indel = start
      
      # special case where the first 9 bases or something of the ref are deleted in the alt, 
      # but the deletion occurs at the first base; this is an in-frame mutation but can't look up
      # the codon position 
      if(start==1){
        LastReferencePosition=0
      }else{
        LastReferencePosition = refDF[(start-1), "codonPos"]
      }
      
      if( ((lengthIndel %% 3) == 0) & (LastReferencePosition %in% c(0,3)) ){
        # check to see if there's a stop codon inserted/deleted
        rows_indel = consistentSNPs %>% filter(POS %in% l)
        rows_indelRef = paste(rows_indel$REF, collapse="")
        rows_indelAlt = paste(rows_indel$ALT, collapse="")
        
        if(grepl("\\*", rows_indelAlt)){
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
        # "Absent" in high if the DELETION is absent from the high group
        PresentOrAbsent_High = "absent"
      }else{
        # "Present" in high if the DELETION is present in the high group
        PresentOrAbsent_High = "present"
      }
  
    rowAnnotated = c(genename, Phenotype, CClabel, patientID, Cluster1, Cluster2,start, vartype, PresentOrAbsent_High)
    rowAnnotatedINDELs = c(genename, Phenotype, CClabel, patientID, Cluster1, Cluster2,start, vartype, PresentOrAbsent_High, lengthIndel,RefBases, AltBases )
    
    AnnotatedSNPs = rbind(AnnotatedSNPs, rowAnnotated)
    Big_Indel_DF=rbind(Big_Indel_DF, rowAnnotatedINDELs)
    colnames(Big_Indel_DF) = c("Gene", "Phenotype", "CC", "Patient","Cluster1", "Cluster2","Position","VariantType", "PresentOrAbsentHigh", "LengthIndel","RefBases","AltBases")
    RemainingSNPs = RemainingSNPs %>% filter(!(POS %in% l))
  }
    
  
  
  # Then, if there's anything left in RemainingSNPs after dealing with the indels
  ###############################################################################
  
  # A SNP position refers to where in the ALIGNMENT a SNP falls; not where in the reference sequence
  # itself the snp falls; need to map the reference positions to their alignment positions
  
  if(nrow(RemainingSNPs) > 0){
    
    refDF = refDF %>% left_join(RemainingSNPs %>% select(POS, REF, ALT), by="POS")
    row.names(refDF) = refDF$POS
    

    for(snpInd in 1:nrow(RemainingSNPs)){
      position = RemainingSNPs[snpInd, "POS"]
      rowinfo = refDF[position,]
      refDF$ALT[refDF$ALT==TRUE] <- "T" 
      codonbases = (refDF %>% filter(AAnum == as.character(rowinfo["AAnum"])))
      referencecodonbases = paste(codonbases$refbase, collapse="")
      variantcodonbases  = paste(if_else(is.na(codonbases$ALT), codonbases$refbase, codonbases$ALT), collapse="")

      old = Biostrings::GENETIC_CODE[[referencecodonbases]]

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
    colnames(AnnotatedSNPs) = c("Gene", "Phenotype", "CC", "Patient", "Cluster1", "Cluster2", "Position", "Type", "PresentOrAbsentHigh")
    AnnotatedSNPs = AnnotatedSNPs %>% filter(Type!="Synonymous")
    
    
  
  }
  if(nrow(AnnotatedSNPs) > 0 ){
    AnnotatedSNPs$NumTotalSeqVariants = nrow(AnnotatedSNPs)
    colnames(AnnotatedSNPs) = c("Gene", "Phenotype", "CC", "Patient", "Cluster1", "Cluster2", "Position", "Type", "PresentOrAbsentHigh", "NumTotalVariants")
    
  }
  BigDF = rbind(BigDF, AnnotatedSNPs)
  colnames(BigDF) = c("Gene", "Phenotype", "CC", "Patient", "Cluster1", "Cluster2", "Position", "Type", "PresentOrAbsentHigh", "NumTotalVariants")
}




# Look at copy number variants 

CNV_Variants = data.frame()
for (filen in list.files("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/CNV_Info/")){
  fpath = paste0("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/CNV_Info/", filen )
  CNVs = read.csv(fpath)
  filenamestring = str_remove_all(filen, "_duplications.csv")
  
  filenamestringList = str_split(filenamestring, "_")[[1]]
  Patient= filenamestringList[2]
  CC= filenamestringList[3]
  Cluster1=filenamestringList[4]
  Cluster2=filenamestringList[5]
  Phenotype=filenamestringList[6]
  Gene = paste(filenamestringList[7:length(filenamestringList)], collapse="_")
  
  
  
  HighCopyNumbers = (CNVs %>% filter(HighOrLow=="High"))$NumCopies
  LowCopyNumbers = (CNVs %>% filter(HighOrLow=="Low"))$NumCopies
  if(length(unique(HighCopyNumbers))==1 & length(unique(LowCopyNumbers))==1  ){
    highrep = HighCopyNumbers[1]
    lowrep = LowCopyNumbers[1]
    if(highrep!=lowrep){
      Type = "HighCopy"
      PresentAbsentHigh = if_else(highrep > lowrep, "present", "absent")
      CNV_Variants = rbind(CNV_Variants, c(Gene, Patient, CC, Cluster1, Cluster2,"HigherCopyNumber",Phenotype, PresentAbsentHigh))
      
      # "Gene", "Phenotype", "CC", "Patient", "Cluster1", "Cluster2", "Position", "Type", "PresentOrAbsentHigh", "NumTotalVariants"

      colnames(CNV_Variants)= c("Gene","Patient", "CC","Cluster1", "Cluster2","Type", "Phenotype", "PresentOrAbsentHigh")
    }
  }
}

write.csv(CNV_Variants, "Documents/DataInputGithub/data/IntraPatient/CNV_Variants_Consistent.csv")

write.csv(BigDF, "Documents/DataInputGithub/data/IntraPatient/SNP_Variants_Consistent.csv")
write.csv(Big_Indel_DF, "Documents/DataInputGithub/data/IntraPatient/INDEL_variants_Consistent.csv")

NonSynon = (BigDF %>% filter(Type!="Synonymous"))

test_table = table(BigDF %>% filter(Type!="Synonymous") %>% select( Gene, Phenotype, Patient))

test_table = table(BigDF %>% filter(Type!="Synonymous") %>% select( Gene, Patient) %>% unique() %>% select(Gene))


NonSynonGenesFixed = sapply(NonSynon$Gene,function(x) str_replace(x, "extdb_", "extdb:"))
NonSynonGenesFixed = sapply(NonSynonGenesFixed,function(x) str_replace(x, "_", "-"))

Annotations = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/RoaryResultsPGAP2022/gene_presence_absence_new_WithPanGenomeIDs.csv") %>% select(Gene, Annotation, Non.unique.Gene.name, Avg.group.size.nuc) 


Annotations$FixedGeneNames = sapply(Annotations$Gene, function(x) str_replace(x,"\\'", "_"))
Annotations$FixedGeneNames = sapply(Annotations$FixedGeneNames, function(x) str_replace(x,"\\(", "_"))
Annotations$FixedGeneNames = sapply(Annotations$FixedGeneNames, function(x) str_replace(x,"\\)", "_"))
Annotations$FixedGeneNames = sapply(Annotations$FixedGeneNames, function(x) str_replace(x,":", "_"))
Annotations$FixedGeneNames = sapply(Annotations$FixedGeneNames, function(x) str_replace(x,"-", "_"))
annotatd = Annotations %>% filter(FixedGeneNames %in% NonSynon$Gene)



BigDF %>% filter(Type!="Synonymous") %>% filter(Gene %in% c("group_2378","coa" ))

BigDF %>% filter(Type!="Synonymous") %>% filter(Gene %in% c("group_70","deoC" ))

BigDF %>% filter(Type!="Synonymous") %>% filter(Gene=="extdb_pgaptmp_002685")

BigDF %>% filter(Type!="Synonymous") %>% filter(Gene %in% c("extdb_pgaptmp_001880","group_237" ))

BigDF %>% filter(Type!="Synonymous") %>% filter(Gene %in% c("group_769","group_2438" ))

multiple_patients = names(test_table[test_table>1])
