#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/ST20.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1442.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1462.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1483.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1490.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1509.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1496.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1471.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1430.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1473.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/ST20.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/RoaryOutput/core_gene_alignment.aln -n 166_CC20.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/List_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/RoaryOutput/core_gene_filelist.txt 10

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/Trees/RAxML_bestTree.166_CC20.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_166_CC20/Trees/patient_166_CC20_clonalframeML.newick
