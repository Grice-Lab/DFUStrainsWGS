#!/bin/bash
# Amy Campbell
# January 2021 
# Map unitigs from dbgwas run to /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_WithAnnotations/Unitigs.fasta
source /home/acampbe/software/miniconda3/bin/activate DBGWASEnv

ReferenceGenome="/home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/SA_502A.fasta"
ReferencePrefix="SA_502A"
UnitigsList="/home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_WithAnnotations/Unitigs.fasta"

#bwa index -p $ReferencePrefix $ReferenceGenome 

bwa mem -O 3 -B 2 $ReferencePrefix $UnitigsList > UnitigAlignment.sam

