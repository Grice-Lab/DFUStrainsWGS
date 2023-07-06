#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC72_CN1.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN691.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN701.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN672.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC72_CN1.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/RoaryOutput/core_gene_alignment.aln -n 131_CC72.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/DFUStrainsList_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/RoaryOutput/core_gene_filelist.txt 4

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/Trees/RAxML_bestTree.131_CC72.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_131_CC72/Trees/patient_131_CC72_clonalframeML.newick
