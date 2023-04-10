# Amy Campbell
# March 2023
# Take simulated reads and map them to genetic markers for CC1 and CC15 in Patient111

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


markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/Patient111GeneMarkers.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/Alignments/"

mkdir -p $outputpath

# Build bowtie index of the marker genes
########################################
bowtie2-build $markergenepath Patient111Markers
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment of simulated reads to 


# Run bowtie2 local alignment on the simulated reads with default params for "very sensitive local"
####################################################################################################

#################
# S. epidermidis
#################
OutputSamFile="SEpi_Markers111.sam"
OutputBamFile="SEpi_Markers111.bam"
OutputBamFileSorted="SEpi_Markers111_sorted.bam"
OutputBamFileSortedBai="SEpi_Markers111_sorted.bam.bai"

bowtie2 --very-sensitive-local -x Patient111Markers -U $S_epi_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_Epidermidis_Coverage_111markers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

#################
# S. pettenkefori
#################

OutputSamFile="SPet_Markers111.sam"
OutputBamFile="SPet_Markers111.bam"
OutputBamFileSorted="SPet_Markers111_sorted.bam"
OutputBamFileSortedBai="SPet_Markers111_sorted.bam.bai"

bowtie2 --very-sensitive-local -x Patient111Markers -U $S_petten_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_Pettenkefori_Coverage_111markers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

###########
# S. aureus
###########

OutputSamFile="Sau_Markers111.sam"
OutputBamFile="Sau_Markers111.bam"
OutputBamFileSorted="Sau_Markers111_sorted.bam"
OutputBamFileSortedBai="Sau_Markers111_sorted.bam.bai"

#bowtie2 --very-sensitive-local -x Patient111Markers -U $S_aureus_readspath -S $outputpath$OutputSamFile

#samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

#samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

#samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_aureus_Coverage_111markers.tsv"

#python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV


