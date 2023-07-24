#!/bin/bash
# Amy Campbell
# Make a big multifasta of phage and plasmid represenative sequences

phagepath=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/ClusteredPhagesRenamed.fasta
plasmidpath=/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/PlasmidReps.fasta
newpath=/home/acampbe/DFU/data/WGS_2020/RoaryResults2022/PhagesPlasmids.fasta

cp $phagepath $newpath
cat $plasmidpath >> $newpath
