#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC5_Mu50.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1447.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1465.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1455.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC5_Mu50.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/RoaryOutput/core_gene_alignment.aln -n 165_CC5.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/DFUStrainsList_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/RoaryOutput/core_gene_filelist.txt 4

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/Trees/RAxML_bestTree.165_CC5.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_165_CC5/Trees/patient_165_CC5_clonalframeML.newick
