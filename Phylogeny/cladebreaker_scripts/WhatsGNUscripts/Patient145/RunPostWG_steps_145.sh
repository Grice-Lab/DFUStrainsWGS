#!bin/bash
# Amy Campbell
# March 2023
# calls different scripts involved in parsing and making use of WhatsGNU output 
# for choosing marker genes for metagenomic shotgun data 

VolcanoPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/Patient145_WhatsGNU_volcano.txt

# Conda environment
source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

# Make list of orthologs to consider based on whatsGNU plotting output
#####################################################################
#Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/WhatsGNUoutput.R $VolcanoPath


OrthoList=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/Patient145_ortholog_list.txt
clusteredprots=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/roary_alignment/results/clustered_proteins
WGrep=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/WG_reports/DORN1151_WhatsGNU_report.txt
Genomefaa=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/DORN1151/annotation/DORN1151.faa

#PanGenomePath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/roary_alignment/results/pan_genome_reference.fa

# Make list of orthologs that aren't within 60% sequence identity to any other orthologous groups
#################################################################################################
#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/UniqueWithinGenome.py $OrthoList $clusteredprots $WGrep $Genomefaa

FFNfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/FFNs

#mkdir -p $FFNfolder

# Move all the necessary .ffns into a folder
############################################
#for f in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/WG_reports/* ; do
	
#	bname=$(basename $f)
#	oldext="_WhatsGNU_report.txt"
#	noext=""
#	genomename=${bname/$oldext/$noext}
#	folder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/$genomename/annotation/
#	cp $folder$genomename".ffn" $FFNfolder/
#done


roaryCSVpath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/roary_alignment/results/gene_presence_absence.csv
UniqueGenes=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/Patient145_unique.txt
otptScript="/home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/"
ReferenceGenome="DORN1151"

# Make multifastas of all the remaining orthologs (if indeed they are found in all the cases and controls by roary) and fastas of one reference genome's version of the gene

#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/OutputMultiFastas_Gene.py $FFNfolder $roaryCSVpath $UniqueGenes $otptScript $ReferenceGenome
#mv ../Patient145_MafftAlignGenes.sh .
#mv ../Patient145_SNPsites.sh .

# Run MAFFT on each of those genes 
##################################
#bsub -e mafft145.e -o mafft145.o -J "patient_145_MAFFT" sh Patient145_MafftAlignGenes.sh

# Run SNP-SITES on each of those genes
######################################
#bsub -e snps145.e -o snps145.o -w "done(patient_145_MAFFT)" sh Patient145_SNPsites.sh

# Extract the variant sites
###########################

refGenome="DORN1151"
caseControl="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/case_controls_patient145.csv"
dirpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/GeneFastas/MAFFT"
otptpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/Patient145_Variants.tsv"

Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/ReadVCFs.R $refGenome $caseControl $dirpath $otptpath

# Make a multifasta of each marker gene from the reference 'case' genome for metagenomic read alignment
#######################################################################################################

GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/GeneFastas"
FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/Patient145_Variants_FinalGeneList.txt"
PatientID="Patient145"
OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/Patient145GeneMarkers.fasta"
ReferenceGenome="DORN1151"
python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

mkdir -p /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient145
cp $OutputfilePath /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient145/
