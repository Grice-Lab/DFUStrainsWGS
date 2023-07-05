# Amy Campbell
# 2023
# First step making intra-patient, intra-CC trees to map phenotype & 
# HGT elements onto

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
  linelistShell=c("#!bin/bash", "# Making intrapatient, intra-CC tree", "mamba ~/mambaforge/bin/activate RoaryEnvNewest")
  
  FolderName=paste0("Patient_",combo )
  outputdirectory=paste0(Outermostfolder, "/",FolderName)
  gffdirectory=paste0(outputdirectory, "/","gffs/")
  
  linelistShell= append(linelistShell, paste0("mkdir ", outputdirectory))
  linelistShell= append(linelistShell, paste0("mkdir ", gffdirectory))
  
  
  
  # paste0("cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/", (SubsetDF$NameReferenceGFF)[1])
}
       

# lineListShell = c("#!bin/bash","# Performing KMC operations for clade-specific markers","", "source ~/mambaforgeOLD/bin/activate ~/mambaforgeOLD/envs/pankmer", "")
# lineListShell = append(lineListShell, paste0("kmc_tools complex ",IntersectOpsFileID ))
# lineListShell = append(lineListShell, paste0("kmc_tools complex ",UnionOpsFileID ))
# lineListShell = append(lineListShell, paste0("kmc_tools simple ",IntersectOutput," ",UnionOutput, " kmers_subtract ",UniqueOutput ))
# lineListShell = append(lineListShell, paste0("rm ",IntersectOutput,"*" ))
# lineListShell = append(lineListShell, paste0("rm ",UnionOutput, "*"))
# lineListShell = append(lineListShell, "")
# lineListShell = append(lineListShell, paste0("kmc_tools transform ", UniqueOutput, " dump ",UniqueOutputText))
# 
# writeLines(lineListShell, ShellScriptOutput)