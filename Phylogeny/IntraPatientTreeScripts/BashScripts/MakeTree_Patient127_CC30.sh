#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC30_MRSA252.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN596.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN613.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN620.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC30_MRSA252.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/RoaryOutput/core_gene_alignment.aln -n 127_CC30.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/List_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/RoaryOutput/core_gene_filelist.txt 4

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/Trees/RAxML_bestTree.127_CC30.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_127_CC30/Trees/patient_127_CC30_clonalframeML.newick
