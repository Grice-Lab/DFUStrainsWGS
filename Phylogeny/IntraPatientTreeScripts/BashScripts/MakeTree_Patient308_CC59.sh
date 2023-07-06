#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC59.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN915.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN968.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN962.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC59.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/RoaryOutput/core_gene_alignment.aln -n 308_CC59.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/List_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/RoaryOutput/core_gene_filelist.txt 4

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/Trees/RAxML_bestTree.308_CC59.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_308_CC59/Trees/patient_308_CC59_clonalframeML.newick -xmfa_file true
