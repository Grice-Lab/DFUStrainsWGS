#!bin/bash
# Amy Campbell
# March 2023
# calls different scripts involved in parsing and making use of WhatsGNU output 
# for choosing marker genes for metagenomic shotgun data 

VolcanoPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/Patient159_WhatsGNU_volcano.txt

# Conda environment
source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

# Make list of orthologs to consider based on whatsGNU plotting output
#####################################################################
#Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/WhatsGNUoutput.R $VolcanoPath


OrthoList=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/Patient159_ortholog_list.txt
clusteredprots=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/roary_alignment/results/clustered_proteins
WGrep=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/WG_reports/DORN1364_WhatsGNU_report.txt
Genomefaa=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/DORN1364/annotation/DORN1364.faa

#PanGenomePath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/roary_alignment/results/pan_genome_reference.fa

# Make list of orthologs that aren't within 60% sequence identity to any other orthologous groups
#################################################################################################
#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/UniqueWithinGenome.py $OrthoList $clusteredprots $WGrep $Genomefaa

FFNfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/FFNs

#mkdir -p $FFNfolder

# Move all the necessary .ffns into a folder
############################################
#for f in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/WG_reports/* ; do
	
#	bname=$(basename $f)
#	oldext="_WhatsGNU_report.txt"
#	noext=""
#	genomename=${bname/$oldext/$noext}
#	folder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/$genomename/annotation/
#	cp $folder$genomename".ffn" $FFNfolder/
#done


roaryCSVpath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/roary_alignment/results/gene_presence_absence.csv
UniqueGenes=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/Patient159_unique.txt
otptScript="/home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/"
# Make multifastas of all the remaining orthologs (if indeed they are found in all the cases and controls by roary) and fastas of one reference genome's version of the gene
ReferenceGenome="DORN1364"

#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/OutputMultiFastas_Gene.py $FFNfolder $roaryCSVpath $UniqueGenes $otptScript $ReferenceGenome
#mv ../Patient159_MafftAlignGenes.sh .
#mv ../Patient159_SNPsites.sh .

# Run MAFFT on each of those genes 
##################################
#bsub -e mafft159.e -o mafft159.o -J "patient_159_MAFFT" sh Patient159_MafftAlignGenes.sh

# Run SNP-SITES on each of those genes
######################################
#bsub -e snps159.e -o snps159.o -w "done(patient_159_MAFFT)" sh Patient159_SNPsites.sh

# Extract the variant sites
###########################

refGenome="DORN1364"
caseControl="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/case_controls_patient159.csv"
dirpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/GeneFastas/MAFFT"
otptpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/Patient159_Variants.tsv"

Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/ReadVCFs.R $refGenome $caseControl $dirpath $otptpath

# Make a multifasta of each marker gene from the reference 'case' genome for metagenomic read alignment
#######################################################################################################

GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/GeneFastas"
FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/Patient159_Variants_FinalGeneList.txt"
PatientID="Patient159"
OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/Patient159GeneMarkers.fasta"
ReferenceGenome="DORN1364"

python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

mkdir -p /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159
cp $OutputfilePath /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159/
