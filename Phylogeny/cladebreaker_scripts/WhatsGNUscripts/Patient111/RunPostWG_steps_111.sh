#!bin/bash
# Amy Campbell
# March 2023
# calls different scripts involved in parsing and making use of WhatsGNU output 
# for choosing marker genes for metagenomic shotgun data 

VolcanoPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/Patient111_WhatsGNU_volcano.txt

# Conda environment
source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

# Make list of orthologs to consider based on whatsGNU plotting output
#####################################################################
#Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/WhatsGNUoutput.R $VolcanoPath


OrthoList=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/Patient111_ortholog_list.txt
clusteredprots=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/roary_alignment/results/clustered_proteins
WGrep=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/WG_reports/DORN152_WhatsGNU_report.txt
Genomefaa=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/DORN152/annotation/DORN152.faa

#PanGenomePath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/roary_alignment/results/pan_genome_reference.fa

# Make list of orthologs that aren't within 60% sequence identity to any other orthologous groups
#################################################################################################
#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/UniqueWithinGenome.py $OrthoList $clusteredprots $WGrep $Genomefaa

FFNfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/FFNs

#mkdir -p $FFNfolder

# Move all the necessary .ffns into a folder
############################################
#for f in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/WG_reports/* ; do
	
#	bname=$(basename $f)
#	oldext="_WhatsGNU_report.txt"
#	noext=""
#	genomename=${bname/$oldext/$noext}
#	folder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/$genomename/annotation/
#	cp $folder$genomename".ffn" $FFNfolder/
#done


roaryCSVpath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/roary_alignment/results/gene_presence_absence.csv
UniqueGenes=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/Patient111_unique.txt
otptScript="/home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/"
ReferenceGenome="DORN152"

# Make multifastas of all the remaining orthologs (if indeed they are found in all the cases and controls by roary) and fastas of one reference genome's version of the gene

#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/OutputMultiFastas_Gene.py $FFNfolder $roaryCSVpath $UniqueGenes $otptScript $ReferenceGenome
#mv ../Patient111_MafftAlignGenes.sh .
#mv ../Patient111_SNPsites.sh .

# Run MAFFT on each of those genes 
##################################
#bsub -e mafft111.e -o mafft111.o -J "patient_111_MAFFT" sh Patient111_MafftAlignGenes.sh

# Run SNP-SITES on each of those genes
######################################
#bsub -e snps111.e -o snps111.o -w "done(patient_111_MAFFT)" sh Patient111_SNPsites.sh

# Extract the variant sites
###########################

refGenome="DORN152"
caseControl="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/case_controls_patient111.csv"
dirpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/GeneFastas/MAFFT"
otptpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/Patient111_Variants.tsv"

Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/ReadVCFs.R $refGenome $caseControl $dirpath $otptpath

# Make a multifasta of each marker gene from the reference 'case' genome for metagenomic read alignment
#######################################################################################################

GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/GeneFastas"
FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/Patient111_Variants_FinalGeneList.txt"
PatientID="Patient111"
OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/Patient111GeneMarkers.fasta"
ReferenceGenome="DORN152"
python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

mkdir -p /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111
cp $OutputfilePath /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/
