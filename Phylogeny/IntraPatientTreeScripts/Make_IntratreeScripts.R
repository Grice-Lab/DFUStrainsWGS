# Amy Campbell
# 2023
# First step making intra-patient, intra-CC trees to map phenotype & 
# HGT elements onto. Make roary core gene alignment, raxml tree

isolateinfo = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/DFU_Staph_aureus_isolates.csv")
isolateinfo$DORN = paste0("DORN", isolateinfo$Doern.lab.bank.)
isolates = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/Phylogeny2022Data/CCMapPlotting.csv")
CCreferencenames = read.csv("/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/IntraPatientRefs.csv")[, 1:2]

# Other lines to add to each bash script
Outermostfolder="/home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/"

isolates = isolates %>% filter(! (DORN %in% c("DORN1176", "DORN429", "DORN685", "DORN946")))

isolates = isolates %>% left_join(isolateinfo %>% select(DORN, patient_id ), by="DORN")
isolates$id_cc = paste0(isolates$patient_id,"_", isolates$CCLabel)

# unique combinations of patient and CC, and how many isolates fall into each
patientCCtable = table(isolates$id_cc)

AtLeast3 = patientCCtable[patientCCtable>=3]

FilteredIsolatesAtLeast3 = isolates %>% filter(id_cc %in% names(AtLeast3))

# For every unique Patient ID / CC combination, 
# Make a bash script that runs roary and raxML on it 
# And then runs ClonalFrameML on the resulting tree
FilteredIsolatesAtLeast3 = FilteredIsolatesAtLeast3 %>% left_join(CCreferencenames, by="CCLabel")
for (combo in unique(FilteredIsolatesAtLeast3$id_cc)){
  SubsetDF = FilteredIsolatesAtLeast3 %>% filter(id_cc==combo)
  ShellScriptOutput=paste0("Phylogeny/IntraPatientTreeScripts/BashScripts/", "MakeTree_Patient", combo, ".sh")
  linelistShell=c("#!bin/bash", "# Making intrapatient, intra-CC tree", "mamba ~/mambaforge/bin/activate RoaryEnvNewest", "")
  linelistShell_CFML = c("#!bin/bash", "# Fixing branch lengths of intrapatient, intra-CC tree using ClonalFrameML", "mamba ~/mambaforge/bin/activate RoaryEnvNewest", "")
  
  
  FolderName=paste0("Patient_",combo )
  outputdirectory=paste0(Outermostfolder, FolderName)
  gffdirectory=paste0(outputdirectory, "/","gffs/")
  Treedirectory=paste0(outputdirectory, "/","Trees/")
  NumGenomes = nrow(SubsetDF) + 1 # add one for the reference genome 
    
  linelistShell= append(linelistShell, paste0("mkdir -p ", outputdirectory))
  linelistShell= append(linelistShell, paste0("mkdir -p ", gffdirectory))
  linelistShell=append(linelistShell, paste0("mkdir -p ", Treedirectory))
  
  linelistShell= append(linelistShell, "")
  linelistShell = append (linelistShell, paste0("cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/", (SubsetDF$NameReferenceGFF)[1], ".gff ", gffdirectory))
  
  for(DORNstring in SubsetDF$DORN){
    linelistShell = append (linelistShell, paste0("cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/", paste0(DORNstring, ".gff ", gffdirectory)))
    
  }
 
  linelistShell = append (linelistShell, paste0("cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/", (SubsetDF$NameReferenceGFF)[1], ".gff ", gffdirectory))
  
  RoaryOutput=paste0(outputdirectory, "/RoaryOutput")
  RoaryCall = paste0("roary -e -z -p -4 -f ", RoaryOutput, " ", gffdirectory, "*")
  
  linelistShell = append(linelistShell, RoaryCall)
  linelistShell = append(linelistShell, "")
  
  coregenealignment=paste0(RoaryOutput, "/core_gene_alignment.aln")
  TreeName=paste0(combo,".newick")
  
  linelistShell=append(linelistShell, paste0("raxmlHPC -m GTRGAMMA -p 19104 -s ", coregenealignment," -n ", TreeName))
  linelistShell=append(linelistShell, paste0("mv *newick* ", Treedirectory))
  
  BestTree = paste0(Treedirectory,"RAxML_bestTree.",TreeName)
  coregenelistpath=paste0(RoaryOutput, "/core_gene_filelist.txt")
  ListCoreGenesCall=paste0("Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/List_Core_Alignment_Files.R ", RoaryOutput, "/ ", coregenelistpath, " ", NumGenomes)
  MakeXMFACAll=paste0("python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py ", RoaryOutput, " ",coregenelistpath )
  
  XMFApath=paste0(RoaryOutput, "/core_genes.xmfa")
  CFMLtreeName=paste0(Treedirectory,"patient_", combo, "_clonalframeML.newick")
  CFML_Call = paste0("ClonalFrameML ", BestTree, " ", XMFApath, " ",CFMLtreeName, " -xmfa_file true")
  
  linelistShell = append(linelistShell, "")
  linelistShell = append(linelistShell, "# Prepare core gene-by-gene alignment for input into ClonalFrameML")
  linelistShell = append(linelistShell, "##################################################################")
  linelistShell = append(linelistShell, ListCoreGenesCall)
  linelistShell = append(linelistShell, "")
  linelistShell = append(linelistShell, MakeXMFACAll)
  linelistShell = append(linelistShell, "")
  linelistShell = append(linelistShell, "# Call ClonalFrameML")
  linelistShell = append(linelistShell, "#####################")
  linelistShell = append(linelistShell, CFML_Call)
  writeLines(linelistShell, ShellScriptOutput)
  
}

