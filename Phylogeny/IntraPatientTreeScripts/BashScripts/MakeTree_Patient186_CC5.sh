#!bin/bash
# Making intrapatient, intra-CC tree
mamba ~/mambaforge/bin/activate RoaryEnvNewest

mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
mkdir -p /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/Trees/

cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC5_Mu50.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1829.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1844.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1849.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1869.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1902.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1952.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN1968.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN2004.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN2034.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN2075.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN2127.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN2139.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN2166.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/DORN2181.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
cp /home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/CC5_Mu50.gff /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/
roary -e -z -p -4 -f /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/gffs/*

raxmlHPC -m GTRGAMMA -p 19104 -s /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/RoaryOutput/core_gene_alignment.aln -n 186_CC5.newick
mv *newick* /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/Trees/

# Prepare core gene-by-gene alignment for input into ClonalFrameML
##################################################################
Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/DFUStrainsList_Core_Alignment_Files.R /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/RoaryOutput/ /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/RoaryOutput/core_gene_filelist.txt 15

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/Create_XMFA_File_Generalized.py /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/RoaryOutput /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/RoaryOutput/core_gene_filelist.txt

# Call ClonalFrameML
#####################
ClonalFrameML /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/Trees/RAxML_bestTree.186_CC5.newick /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/RoaryOutput/core_genes.xmfa /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/Patient_186_CC5/Trees/patient_186_CC5_clonalframeML.newick
