# Amy Campbell
# March 2023
# Take simulated reads and map them to genetic markers for CC1 and CC15 in Patient124

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


markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/Patient124GeneMarkers.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/Alignments/"

mkdir -p $outputpath

# Build bowtie index of the marker genes
########################################
bowtie2-build $markergenepath Patient124Markers
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment of simulated reads to 


# Run bowtie2 local alignment on the simulated reads with default params for "very sensitive local"
####################################################################################################

#################
# S. epidermidis
#################
OutputSamFile="SEpi_Markers124.sam"
OutputBamFile="SEpi_Markers124.bam"
OutputBamFileSorted="SEpi_Markers124_sorted.bam"
OutputBamFileSortedBai="SEpi_Markers124_sorted.bam.bai"

bowtie2 --very-sensitive-local -x Patient124Markers -U $S_epi_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_Epidermidis_Coverage_124markers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

rm $outputpath$OutputBamFileSortedBai
rm $outputpath$OutputBamFileSorted
rm $outputpath$OutputBamFile
rm $outputpath$OutputSamFile

#################
# S. pettenkefori
#################

OutputSamFile="SPet_Markers124.sam"
OutputBamFile="SPet_Markers124.bam"
OutputBamFileSorted="SPet_Markers124_sorted.bam"
OutputBamFileSortedBai="SPet_Markers124_sorted.bam.bai"

bowtie2 --very-sensitive-local -x Patient124Markers -U $S_petten_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_Pettenkefori_Coverage_124markers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV
rm $outputpath$OutputBamFileSortedBai
rm $outputpath$OutputBamFileSorted
rm $outputpath$OutputBamFile
rm $outputpath$OutputSamFile
###########
# S. aureus
###########

OutputSamFile="Sau_Markers124.sam"
OutputBamFile="Sau_Markers124.bam"
OutputBamFileSorted="Sau_Markers124_sorted.bam"
OutputBamFileSortedBai="Sau_Markers124_sorted.bam.bai"

#bowtie2 --very-sensitive-local -x Patient124Markers -U $S_aureus_readspath -S $outputpath$OutputSamFile

#samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

#samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

#samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_aureus_Coverage_124markers.tsv"

#python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV


