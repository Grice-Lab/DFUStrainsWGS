# Amy Campbell
# March 2023

# Similar to ReadVCFs.R but specifically designed for Patient 124 who has 3 strains 
# Takes in:
# 1) Path to folder containing the VCFs 
# 2) case control csv which was input into whatsGNU
# 3) Name of the reference genome

library(stringr)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)

# takes in: 
##############
ReferenceGenome=args[1]

DirectoryPath=args[2]

OutputPath= args[3]



CC5List = c("DORN781","GCA_000555645.1", "GCA_000530705.1", "GCA_000554505.1", "GCA_002982825.1")
CC9List=c("DORN653", "DORN659", "GCA_900081385.1","GCA_003813255.1","GCA_900082315.1")
CC8List = c("DORN717", "DORN657", "GCA_001353315.1", "GCA_000541135.1", "GCA_000017085.1")
AllGenomes=c(CC5List,CC9List, CC8List )

FilesInput = list.files(DirectoryPath)
FilesInput = FilesInput[grepl("SNPs_",FilesInput )]


FullVariantDF=data.frame()
FinalGeneList =c()


for(f in FilesInput){

  fullpath=file.path(DirectoryPath, f)
  
  # name of gene
  InputOrthologName=str_replace(f, "SNPs_", "")
  
  # Read in the VCF
  inputDF = read.table(fullpath,skip=3,sep='\t', comment.char="$",header=T)
  
  # Fixing stupid read table assumption that string T=true lol 
  inputDF = inputDF %>% mutate(REF = if_else(REF==TRUE, "T", as.character(REF)))
  inputDF = inputDF %>% mutate(ALT = if_else(ALT==TRUE, "T", as.character(ALT)))
  
  # Remove the extraneous instance of the reference gene which we used for aligning all of them 
  inputDF = inputDF %>% select(-paste0(ReferenceGenome, ".1"))
  
  # Get the length of the alignment
  inputstring=file(fullpath, open="r")
  inputstringLength=(readLines(inputstring))[2]
  lengthAlign = (str_split(inputstringLength,"length=" ))[[1]][2]
  lengthAlign = as.integer(str_split(lengthAlign, ">")[[1]][1])
  close(inputstring)
  
  
  # If >10% of the length of the whole alignment is represented in variants, OR if there's columns in the SNPs file that aren't in cases or controls
  if( (nrow(inputDF) > .1*lengthAlign) | length(setdiff(AllGenomes, colnames(inputDF))) >0 ){
    print(InputOrthologName)
    print("Too many mismatches or SNP sites couldn't be generated for all genomes. Remove this gene :(")
  }else{
    # Otherwise, proceed. 
    inputDF = inputDF %>% filter(nchar(ALT)==1)
    
    # If there are more than 0 variants that are biallelic, 
    if(nrow(inputDF) > 0){
      
      # Indices for site type A (a variant is present in all of CC9 genomes but neither CC8 nor CC5)
      CC5Cols = inputDF %>% select(CC5List)
      CC9Cols= inputDF %>% select(CC9List)
      CC8Cols= inputDF %>% select(CC8List)
      
      # Variants found in none of the CC5 cols 
      NoCC5 = which(rowSums(CC5Cols)==0)
      # Variants found in none of the CC9 cols 
      NoCC9 = which(rowSums(CC9Cols)==0)
      # Variants found in none of the CC8 cols 
      NoCC8 = which(rowSums(CC8Cols)==0)
      # Variants found in all of the CC9 cols 
      AllCC9= which(rowSums(CC9Cols)==ncol(CC9Cols))
      # Variants found in all of the CC9 cols
      AllCC8 = which(rowSums(CC8Cols)==ncol(CC8Cols))
      

      # In 'case A,' variants which are present in all CC9 but no CC8 and no CC5
      CaseA_indices = intersect(intersect(NoCC8, NoCC5), AllCC9)
      
      # In 'case B,' variants which are present in all CC8 but no CC9 or CC5
      CaseB_indices = intersect(intersect(NoCC9, NoCC5), AllCC8)
      
      
      RowsUseA = inputDF[CaseA_indices, ] %>% select(REF, ALT, POS)
      RowsUseB = inputDF[CaseB_indices, ] %>% select(REF, ALT, POS)
      
      # If there's variants for both cases, 
      if( (length(CaseA_indices)>0) & (length(CaseB_indices)>0)) {
        RowsUseA$Gene = InputOrthologName
        colnames(RowsUseA) = c("Ref", "Alt", "Position", "Gene")
        RowsUseA$Type="A"
        
        RowsUseB$Gene = InputOrthologName
        colnames(RowsUseB) = c("Ref", "Alt", "Position", "Gene")
        RowsUseB$Type="B"
        
        AddedDF=rbind(RowsUseA, RowsUseB)
        
        # If there's variants for A only, 
      }else if((length(CaseA_indices)>0) & (length(CaseB_indices)==0)){
        RowsUseA$Gene = InputOrthologName
        colnames(RowsUseA) = c("Ref", "Alt", "Position", "Gene")
        RowsUseA$Type="A"
        AddedDF=RowsUseA
        print(paste0("No variants for case B exist in ", InputOrthologName))
        
        # If there's variants for B only, 
      }else if((length(CaseB_indices)>0) & (length(CaseA_indices)==0)){
        RowsUseB$Gene = InputOrthologName
        colnames(RowsUseB) = c("Ref", "Alt", "Position", "Gene")
        RowsUseB$Type="B"
        AddedDF=RowsUseB
        print(paste0("No variants for case A exist in ", InputOrthologName))
        #  If there's variants for neither, empty dataframe
      }else if( (length(CaseA_indices)==0) & (length(CaseB_indices)==0)){
        AddedDF=data.frame()
        print(paste0("No variants exist for either case in ", InputOrthologName))
      }
      
      FullVariantDF = rbind(FullVariantDF, AddedDF)
    
# Close if statement for at least 1 variant included
    }
    
# Close else statement for there being all the necessary genomes included 
      }
  
# Close for loop
  }

OutputGeneList = str_replace(OutputPath, ".tsv", "_FinalGeneList.txt")

colnames(FullVariantDF) = c("Case", "Control", "Position", "Gene", "Type")
write.table(FullVariantDF, file=OutputPath, quote=F)

write.table(FullVariantDF, file=OutputPath, quote=F,col.names = F, row.names=F)

