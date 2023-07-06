#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC8_Newman.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1006.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN731.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN788.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN848.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN936.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN973.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN892.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC8_Newman.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/RoaryOutput/core_gene_alignment.aln -n 132_CC8.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/DFUStrainsList_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/RoaryOutput/core_gene_filelist.txt 8

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/Trees/RAxML_bestTree.132_CC8.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_132_CC8/Trees/patient_132_CC8_clonalframeML.newick
