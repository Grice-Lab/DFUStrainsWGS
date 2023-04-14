
library(ggtree)
library(ape)
library(phytools)
library(stringr)

# Polyclonal Patients
######################

TreeFolder="~/Documents/DataInputGithub/data/TreesDFUPatients/"
TreeFilesList = list.files(TreeFolder)
Patient_Genome_Info = read.csv("~/Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")
Patient_Genome_Info$DORN = paste0("DORN", Patient_Genome_Info$Doern.lab.bank.)
Patient_Genome_Info_Visit = Patient_Genome_Info %>% select(DORN, patient_id, visit)

Patient_Genome_Info = Patient_Genome_Info %>% select(DORN, patient_id)
Treeplotlist = c()
for(f in 1:length(TreeFilesList)){
  TreeString = str_replace(TreeFilesList[f], "RAxML_bestTree.patient_", "")
  TreeString = str_replace(TreeString, ".newick", "")
  
  Treefile=paste0(TreeFolder, TreeFilesList[f])
  TreeInput= ggtree::read.tree(Treefile)
  RootedTree = phytools::midpoint.root(TreeInput)
  TreePlot = ggtree::ggtree(RootedTree, size=.25, layout = "circular") + geom_tiplab(size=3)
  ggsave(TreePlot, file=paste0("~/Documents/Saureus_Genomics_Paper/CladebreakerPlots/patient_", TreeString, "_tree.pdf"))
}


CCmapping=read.csv("data/Phylogeny2022Data/CCMapPlotting.csv")
CCsPatients = Patient_Genome_Info_Visit %>% left_join(CCmapping, by="DORN")
