#!bin/bash
# Amy Campbell
# March 2023
# calls different scripts involved in parsing and making use of WhatsGNU output 
# for choosing marker genes for metagenomic shotgun data 

# using SaureusCC8_WhatsGNU_report.txt as the case

VolcanoPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/PatientSimulated_WhatsGNU_volcano.txt

# Conda environment
source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGnuEnv2

# Make list of orthologs to consider based on whatsGNU plotting output
#####################################################################
#Rscript ../WhatsGNUoutput.R $VolcanoPath


OrthoList=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/PatientSimulated_ortholog_list.txt
clusteredprots=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/roary_alignment/results/clustered_proteins
WGrep=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/WG_reports/SaureusCC8_WhatsGNU_report.txt
Genomefaa=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/SaureusCC8/annotation/SaureusCC8.faa



FFNfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/FFNs/
#mkdir -p $FFNfolder
#find /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/ -name '*.ffn' -exec cp -prv '{}' $FFNfolder ';'
#rm /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/FFNs/GCA_002310435.1.ffn
#rm /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/FFNs/GCA_000009645.1.ffn

# Make list of orthologs that aren't within 60% sequence identity to any other orthologous groups
#################################################################################################
#python3 ../UniqueWithinGenome.py $OrthoList $clusteredprots $WGrep $Genomefaa

roaryCSVpath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/roary_alignment/results/gene_presence_absence.csv
UniqueGenes=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/PatientSimulated_unique.txt
otptScript="/home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/"
# Make multifastas of all the remaining orthologs (if indeed they are found in all the cases and controls by roary) and fastas of one reference genome's version of the gene

#python3 OutputMultiFastas_Gene_Simulated.py $FFNfolder $roaryCSVpath $UniqueGenes $otptScript


# Run MAFFT on each of those genes 
##################################
#bsub -e mafft.e -o mafft.o -J "patient_simulated_MAFFT" sh PatientSimulated_MafftAlignGenes.sh

# Run SNP-SITES on each of those genes
######################################
#bsub -e snps.e -o snps.o sh PatientSimulated_SNPsites.sh


# Extract the variant sites
###########################

refGenome="SaureusCC8"
caseControl="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/case_controls_simulated.csv"
dirpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/GeneFastas/MAFFT"
otptpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/markergenesSimulatedCC5CC8.tsv"

# List mafft files where ref contains gaps (don't use these)
############################################################
#echo "Missing" > gappedSimulated.txt

#for fname in $dirpath/Mafft_*; do

#        awk 'FNR==2{if(/-/) print FILENAME}' $fname >> gappedSimulated.txt

#done

removePath=$(pwd)/gappedSimulated.txt

#echo $removePath

#Rscript ../ReadVCFs.R $refGenome $caseControl $dirpath $otptpath $removePath

# Make a multifasta of each marker gene from the reference 'case' genome for metagenomic read alignment
#######################################################################################################

GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/GeneFastas"
FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/markergenesSimulatedCC5CC8_FinalGeneList.txt"
PatientID="PatientSimulated"
OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/PatientSimulated/PatientSimulatedGeneMarkers.fasta"
ReferenceGenome="SaureusCC8"
python3 ../MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

