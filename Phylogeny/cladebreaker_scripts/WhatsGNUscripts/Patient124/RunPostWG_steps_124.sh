#!bin/bash
# Amy Campbell
# March 2023
# calls different scripts involved in parsing and making use of WhatsGNU output 
# for choosing marker genes for metagenomic shotgun data 

VolcanoPathCC9=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/Patient124_CC9_WhatsGNU_volcano.txt
VolcanoPathCC8=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC8/Patient124_CC8_WhatsGNU_volcano.txt

# Conda environment
source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

# Make list of orthologs to consider based on whatsGNU plotting output
#####################################################################
#Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/WhatsGNUoutput.R $VolcanoPathCC9
#Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/WhatsGNUoutput.R $VolcanoPathCC8

OrthoListCC9=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/Patient124_CC9_ortholog_list.txt
OrthoListCC8=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC8/Patient124_CC8_ortholog_list.txt
clusteredprots=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/roary_alignment/results/clustered_proteins

WGrep=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/WG_reports/DORN781_WhatsGNU_report.txt
Genomefaa=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/DORN781/annotation/DORN781.faa

# Make list of orthologs that aren't within 60% sequence identity to any other orthologous groups
#################################################################################################
#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/UniqueWithinGenome.py $OrthoListCC9 $clusteredprots $WGrep $Genomefaa
#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/UniqueWithinGenome.py $OrthoListCC8 $clusteredprots $WGrep $Genomefaa

FFNfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/FFNs

#mkdir -p $FFNfolder

# Move all the necessary .ffns into a folder
############################################
#for f in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/WG_reports/* ; do
	
#	bname=$(basename $f)
#	oldext="_WhatsGNU_report.txt"
#	noext=""
#	genomename=${bname/$oldext/$noext}
#	folder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/$genomename/annotation/
#	cp $folder$genomename".ffn" $FFNfolder/
#done




roaryCSVpath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/roary_alignment/results/gene_presence_absence.csv
UniqueGenescc9=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/Patient124_CC9_unique.txt
UniqueGenescc8=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC8/Patient124_CC8_unique.txt

otptScript="/home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/"
ReferenceGenome="DORN781"

# Make multifastas of all the remaining orthologs (if indeed they are found in all the cases and controls by roary) and fastas of one reference genome's version of the gene

#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/OutputMultiFastas_Gene_124.py $FFNfolder $roaryCSVpath $UniqueGenescc9 $UniqueGenescc8 $otptScript $ReferenceGenome
#mv ../Patient124_CC9_MafftAlignGenes.sh .
#mv ../Patient124_CC9_SNPsites.sh .

# Run MAFFT on each of those genes 
##################################
#bsub -e mafft124.e -o mafft124.o -J "patient_124_MAFFT" sh Patient124_CC9_MafftAlignGenes.sh

# Run SNP-SITES on each of those genes
######################################
#bsub -e snps124.e -o snps124.o -w "done(patient_124_MAFFT)" sh Patient124_CC9_SNPsites.sh

# Extract the variant sites
###########################

refGenome="DORN781"
caseControl="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/case_controls_patient124_cc9.csv"
dirpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/GeneFastas/MAFFT"
otptpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/Patient124_Variants.tsv"

Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/ReadVCFs_Patient124.R $refGenome $dirpath $otptpath

# Make a multifasta of each marker gene from the reference 'case' genome for metagenomic read alignment
#######################################################################################################

GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/GeneFastas"
FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/Patient124_Variants_FinalGeneList.txt"
PatientID="Patient124_CC9"
OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9/Patient124GeneMarkers.fasta"
ReferenceGenome="DORN781"
python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

mkdir -p /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124
cp $OutputfilePath /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/

