#!bin/bash
# Amy Campbell
# March 2023
# calls different scripts involved in parsing and making use of WhatsGNU output 
# for choosing marker genes for metagenomic shotgun data 

VolcanoPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/Patient191_WhatsGNU_volcano.txt

# Conda environment
source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

# Make list of orthologs to consider based on whatsGNU plotting output
#####################################################################
#Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/WhatsGNUoutput.R $VolcanoPath


OrthoList=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/Patient191_ortholog_list.txt
clusteredprots=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/roary_alignment/results/clustered_proteins
WGrep=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/WG_reports/DORN1943_WhatsGNU_report.txt
Genomefaa=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/DORN1943/annotation/DORN1943.faa

#PanGenomePath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/roary_alignment/results/pan_genome_reference.fa

# Make list of orthologs that aren't within 60% sequence identity to any other orthologous groups
#################################################################################################
#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/UniqueWithinGenome.py $OrthoList $clusteredprots $WGrep $Genomefaa

FFNfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/FFNs

mkdir -p $FFNfolder

# Move all the necessary .ffns into a folder
############################################
#for f in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/WG_reports/* ; do
	
#	bname=$(basename $f)
#	oldext="_WhatsGNU_report.txt"
#	noext=""
#	genomename=${bname/$oldext/$noext}
#	folder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/$genomename/annotation/
#	cp $folder$genomename".ffn" $FFNfolder/
#done


roaryCSVpath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/roary_alignment/results/gene_presence_absence.csv
UniqueGenes=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/Patient191_unique.txt
otptScript="/home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/"
# Make multifastas of all the remaining orthologs (if indeed they are found in all the cases and controls by roary) and fastas of one reference genome's version of the gene
ReferenceGenome="DORN1943"

#python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/OutputMultiFastas_Gene.py $FFNfolder $roaryCSVpath $UniqueGenes $otptScript $ReferenceGenome
#mv ../Patient191_MafftAlignGenes.sh .
#mv ../Patient191_SNPsites.sh .

# Run MAFFT on each of those genes 
##################################
#bsub -e mafft191.e -o mafft191.o -J "patient_191_MAFFT" sh Patient191_MafftAlignGenes.sh

# Run SNP-SITES on each of those genes
######################################
#bsub -e snps191.e -o snps191.o -w "done(patient_191_MAFFT)" sh Patient191_SNPsites.sh

# Extract the variant sites
###########################

refGenome="DORN1943"
caseControl="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/case_controls_patient191.csv"
dirpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/GeneFastas/MAFFT"
otptpath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/Patient191_Variants.tsv"

Rscript /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/ReadVCFs.R $refGenome $caseControl $dirpath $otptpath

# Make a multifasta of each marker gene from the reference 'case' genome for metagenomic read alignment
#######################################################################################################

GeneFastasDirectory="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/GeneFastas"
FinallGeneListPath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/Patient191_Variants_FinalGeneList.txt"
PatientID="Patient191"
OutputfilePath="/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/Patient191GeneMarkers.fasta"
ReferenceGenome="DORN1943"
python3 /home/acampbe/DFUStrainsWGS/Phylogeny/cladebreaker_scripts/WhatsGNUscripts/MakeCaseReference.py $GeneFastasDirectory $FinallGeneListPath $PatientID $OutputfilePath $ReferenceGenome

mkdir -p /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191
cp $OutputfilePath /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/
