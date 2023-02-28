# Make input files & a list of commands 

library(dplyr)
library(stringr)


ListFiles = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/listStaphFiles.csv", col.names=c("Name"))
Patient_Genome_Info = read.csv("~/Desktop/GriceLabGit/DFUStrainsWGS/data/DFU_Staph_aureus_isolates.csv")
Patient_Genome_Info$Doern.lab.bank. = sapply(Patient_Genome_Info$Doern.lab.bank., as.character)
ListFiles$DOERNno = sapply(ListFiles$Name, function(x) str_replace(x, "_Final.fna", ""))
ListFiles$Doern.lab.bank. = sapply(ListFiles$DOERNno, function(x) str_replace(x, "DORN", ""))


ListFiles = ListFiles %>% left_join(Patient_Genome_Info %>% select(Doern.lab.bank., patient_id), by="Doern.lab.bank.")

PerPatient = ListFiles %>% count(patient_id) %>% filter(n>1)
fullpath="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/InputCSVs/patient_"
for (p in PerPatient$patient_id){
  subset = ListFiles %>% filter(patient_id==p)
  OutputCSV = data.frame(sample=subset$DOERNno, type="assembly", file_1=paste0("/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolatesFNA/", subset$Name), file_2="")
  Outputname = paste0(fullpath, p, "_input.csv")
  write.csv(OutputCSV, quote = F, row.names = F,file=Outputname)
}
