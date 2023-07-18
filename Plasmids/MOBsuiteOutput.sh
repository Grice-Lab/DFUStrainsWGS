# Amy Campbell
# July 2023
# Call MOBoutput.R to aggregate the plasmid detection output
# make a presence/absence CSV
source ~/mambaforge/bin/activate PlasmidEnv

outputpath=/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/Plasmid_Presence_Absence.csv
mobtyperesults=/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/Results
fullgenomesdir=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates
fulldfpath=/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/FullSummaryPlasmids.csv

Rscript MOBoutput.R $mobtyperesults $fullgenomesdir $outputpath $fulldfpath

