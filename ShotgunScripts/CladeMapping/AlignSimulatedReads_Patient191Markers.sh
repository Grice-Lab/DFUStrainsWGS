# Amy Campbell
# March 2023
# Take simulated reads and map them to genetic markers for CC1 and CC15 in Patient191

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv
#~/mambaforge/condabin/conda

# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/BTdb/"
mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

S_epi_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/S_epidermidissorted.fastq"
S_petten_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/S_pettenkoferisorted.fastq"
S_aureus_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/DORN925sorted.fastq"


markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Patient191GeneMarkers.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Alignments/"

mkdir -p $outputpath

# Build bowtie index of the marker genes
########################################
bowtie2-build $markergenepath Patient191Markers
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment of simulated reads to 


# Run bowtie2 local alignment on the simulated reads with default params for "very sensitive local"
####################################################################################################

#################
# S. epidermidis
#################
OutputSamFile="SEpi_Markers191.sam"
OutputBamFile="SEpi_Markers191.bam"
OutputBamFileSorted="SEpi_Markers191_sorted.bam"
OutputBamFileSortedBai="SEpi_Markers191_sorted.bam.bai"

bowtie2 --very-sensitive-local -x Patient191Markers -U $S_epi_readspath -S $outputpath$OutputSamFile
AlignSimulatedReads_Patient191Markers.sh
samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_Epidermidis_Coverage_191markers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

#################
# S. pettenkefori
#################

OutputSamFile="SPet_Markers191.sam"
OutputBamFile="SPet_Markers191.bam"
OutputBamFileSorted="SPet_Markers191_sorted.bam"
OutputBamFileSortedBai="SPet_Markers191_sorted.bam.bai"

bowtie2 --very-sensitive-local -x Patient191Markers -U $S_petten_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_Pettenkefori_Coverage_191markers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

###########
# S. aureus
###########

OutputSamFile="Sau_Markers191.sam"
OutputBamFile="Sau_Markers191.bam"
OutputBamFileSorted="Sau_Markers191_sorted.bam"
OutputBamFileSortedBai="Sau_Markers191_sorted.bam.bai"

#bowtie2 --very-sensitive-local -x Patient191Markers -U $S_aureus_readspath -S $outputpath$OutputSamFile

#samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

#samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

#samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

#OutputTSV="S_aureus_Coverage_191markers.tsv"

#python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV


rm /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Alignments/*.bam
rm /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Alignments/*.sam 
