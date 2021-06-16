# Amy Campbell
# Make a list of filenames (ex: group_3909.fa.aln) in teh form <gene_name>.fa.aln 
# for input to shell script that loops through files, cats their contents into one .xmfa file separated by '=' 
# This will be input into clonalframeML to correct for recombination on the tree 
# first positional argument is 
# Example call: Rscript List_Core_alignment_Files.R /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/ core_gene_filelist.txt 231

library(dplyr)
library(stringr)
# setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0 | length(args) == 1 | length(args) == 2){
  print("Please provide the folder path containing gene_presence_absence.csv as first argument, the output file name (.txt) as the second argument, and the number of isolates included as the third")
  quit()
} else {
  genelistpath=args[1]
  outputname=args[2]
  num_isolates=as.numeric(as.character(args[3]))
}


gene_list = read.csv(paste0(genelistpath, "gene_presence_absence.csv"))

gene_list = gene_list %>% select(Gene, No..isolates)
core_list = gene_list %>% filter(No..isolates == num_isolates)

core_list$Filename = paste0(core_list$Gene, ".fa.aln")
core_list = core_list %>% mutate(Filename = str_replace(Filename, '\'', '_'))
core_list = core_list %>% mutate(Filename = str_replace(Filename, ' ', '_'))
write.table(file=outputname, core_list$Filename, quote = F, row.names=F, col.names=F)


