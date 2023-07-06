#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC1_MW2.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1531.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1562.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1623.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1540.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1546.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC1_MW2.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/RoaryOutput/core_gene_alignment.aln -n 171_CC1.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/List_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/RoaryOutput/core_gene_filelist.txt 6

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/Trees/RAxML_bestTree.171_CC1.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_171_CC1/Trees/patient_171_CC1_clonalframeML.newick -xmfa_file true
