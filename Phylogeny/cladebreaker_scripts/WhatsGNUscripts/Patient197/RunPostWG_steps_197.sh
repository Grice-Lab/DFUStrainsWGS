#!bin/bash
# Amy Campbell
# March 2023
# calls different scripts involved in parsing and making use of WhatsGNU output 
# for choosing marker genes for metagenomic shotgun data 

VolcanoPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/Patient197_WhatsGNU_volcano.txt

# Conda environment
source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

# Make list of orthologs to consider based on whatsGNU plotting output
#####################################################################
#Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/WhatsGNUoutput_LowBurden.R $VolcanoPath


OrthoList=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/Patient197_ortholog_list.txt
clusteredprots=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/roary_alignment/results/clustered_proteins
WGrep=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/WG_reports/DORN2205_WhatsGNU_report.txt
Genomefaa=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/DORN2205/annotation/DORN2205.faa

#PanGenomePath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/roary_alignment/results/pan_genome_reference.fa

# Make list of orthologs that aren't within 60% sequence identity to any other orthologous groups
#################################################################################################
#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/UniqueWithinGenome.py $OrthoList $clusteredprots $WGrep $Genomefaa

FFNfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/FFNs

mkdir -p $FFNfolder

# Move all the necessary .ffns into a folder
############################################
#for f in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/WG_reports/* ; do
	
#	bname=$(basename $f)
#	oldext="_WhatsGNU_report.txt"
#	noext=""
#	genomename=${bname/$oldext/$noext}
#	folder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/$genomename/annotation/
#	cp $folder$genomename".ffn" $FFNfolder/
#done


roaryCSVpath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/roary_alignment/results/gene_presence_absence.csv
UniqueGenes=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/Patient197_unique.txt
otptScript="/home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/"
# Make multifastas of all the remaining orthologs (if indeed they are found in all the cases and controls by roary) and fastas of one reference genome's version of the gene

#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/OutputMultiFastas_Gene.py $FFNfolder $roaryCSVpath $UniqueGenes $otptScript
#mv ../Patient197_MafftAlignGenes.sh .
#mv ../Patient197_SNPsites.sh .

# Run MAFFT on each of those genes 
##################################
#bsub -e mafft197.e -o mafft197.o -J "patient_197_MAFFT" sh Patient197_MafftAlignGenes.sh

# Run SNP-SITES on each of those genes
######################################
#bsub -e snps197.e -o snps197.o -w "done(patient_197_MAFFT)" sh Patient197_SNPsites.sh

# Extract the variant sites
###########################

refGenome="DORN2205"
caseControl="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/case_controls_patient197.csv"
dirpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/GeneFastas/MAFFT"
otptpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/Patient197_Variants.tsv"

Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/ReadVCFs.R $refGenome $caseControl $dirpath $otptpath

# Make a multifasta of each marker gene from the reference 'case' genome for metagenomic read alignment
#######################################################################################################

GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/GeneFastas"
FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/Patient197_Variants_FinalGeneList.txt"
PatientID="Patient197"
OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/Patient197GeneMarkers.fasta"
ReferenceGenome="DORN2205"
python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

mkdir -p /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient197
cp $OutputfilePath /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient197/
