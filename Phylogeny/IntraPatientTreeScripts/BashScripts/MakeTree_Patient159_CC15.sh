#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC15.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1352.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1364.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1384.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1340.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1334.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1339.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC15.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/RoaryOutput/core_gene_alignment.aln -n 159_CC15.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/List_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/RoaryOutput/core_gene_filelist.txt 7

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/Trees/RAxML_bestTree.159_CC15.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_159_CC15/Trees/patient_159_CC15_clonalframeML.newick
