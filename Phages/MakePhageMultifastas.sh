#!bin/bash
# Amy Campbell
# Making multifastas where each entry (or combination of entries) is a diff 
# genome's instance of a phage

source ~/mambaforge/bin/activate DBSCANSWA

inputfolder=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/IndividualFastas_Hclust
outputfolder=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/Hclust_Multifastas
python3 Make_MultifastaPer_PhageCluster.py $inputfolder $outputfolder
