#!/bin/bash
# Amy Campbell
# January 2021 
# Map unitigs from dbgwas run to /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_WithAnnotations/Unitigs.fasta
source /home/acampbe/software/miniconda3/bin/activate DBGWASEnv

# Path variables -- original run (219 isolates)
###############################################
ReferenceGenome="/home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/SA_502A.fasta"
ReferencePrefix="SA_502A"
UnitigsList="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_WithAnnotations/Unitigs.fasta"

# Index the reference genome
#############################
#bwa index -p $ReferencePrefix $ReferenceGenome 

# Align Unitigs against the reference index
###########################################
#bwa mem -O 3 -B 2 $ReferencePrefix $UnitigsList > UnitigAlignment.sam

# Convert the .sam alignment to a .bam file for efficiency
##########################################################
#samtools view -S -b UnitigAlignment.sam > UnitigAlignment.bam



# Path variables -- new run (full 221 isolates, still k=31)
####################################################
ReferencePrefix="SA_502A"
UnitigsList="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_31/Unitigs.fasta"

# Index the reference genome
#############################
# bwa index -p $ReferencePrefix $ReferenceGenome

# Align Unitigs against the reference index
###########################################
bwa mem -O 3 -B 2 $ReferencePrefix $UnitigsList > /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_31/UnitigAlignment.sam

# Convert the .sam alignment to a .bam file for efficiency
##########################################################
samtools view -S -b /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_31/UnitigAlignment.sam > /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_31/UnitigAlignment.bam


