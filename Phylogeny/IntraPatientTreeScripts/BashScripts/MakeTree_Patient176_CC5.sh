#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC5_Mu50.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1643.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1644.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1663.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1695.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1776.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1808.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1818.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1819.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1834.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1858.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1863.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1885.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1923.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1731.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1646.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1602.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1645.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1732.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC5_Mu50.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/RoaryOutput/core_gene_alignment.aln -n 176_CC5.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/DFUStrainsList_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/RoaryOutput/core_gene_filelist.txt 19

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/Trees/RAxML_bestTree.176_CC5.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_176_CC5/Trees/patient_176_CC5_clonalframeML.newick
