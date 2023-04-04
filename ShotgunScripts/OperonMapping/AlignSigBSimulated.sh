source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv

# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomes/StressOperonMarkers/BTdb/"
mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

S_epi_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/S_epidermidissorted.fastq"
S_petten_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/S_pettenkoferisorted.fastq"
S_aureus_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/DORN925sorted.fastq"
S_lug_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/S_lugdunensissorted.fastq"
S_cap_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/S_capitissorted.fastq"
S_haem_readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/sorted/S_haemolyticussorted.fastq"

markergenepath="/home/acampbe/DFU/data/StressOperonMarkers/SigB_LAC.fasta"

outputpath="/home/acampbe/DFU/data/StressOperonMarkers/Alignments/"

mkdir -p $outputpath

# Build bowtie index of the marker genes
########################################
bowtie2-build $markergenepath SigBMarkers
mv *.bt2 $bowtiepath

#################
# S. epidermidis 
#################

OutputSamFile="SEpi_MarkersSigB.sam"
OutputBamFile="SEpi_MarkersSigB.bam"
OutputBamFileSorted="SEpi_MarkersSigB_sorted.bam"
OutputBamFileSortedBai="SEpi_MarkersSigB_sorted.bam.bai"

#bowtie2 --very-sensitive-local -x SigBMarkers -U $S_epi_readspath -S $outputpath$OutputSamFile

#samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

#samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

#samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

#OutputTSV="S_Epidermidis_Coverage_SigBmarkers.tsv"

#python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

#################
# S. pettenkefori
#################

OutputSamFile="SPet_MarkersSigB.sam"
OutputBamFile="SPet_MarkersSigB.bam"
OutputBamFileSorted="SPet_MarkersSigB_sorted.bam"
OutputBamFileSortedBai="SPet_MarkersSigB_sorted.bam.bai"

#bowtie2 --very-sensitive-local -x SigBMarkers -U $S_petten_readspath -S $outputpath$OutputSamFile

#samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

#samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

#samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

#OutputTSV="S_Pettenkefori_Coverage_SigBmarkers.tsv"

#python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

# S aureus CC1
#########
OutputSamFile="Sau_MarkersSigB.sam"
OutputBamFile="Sau_MarkersSigB.bam"
OutputBamFileSorted="Sau_MarkersSigB_sorted.bam"
OutputBamFileSortedBai="Sau_MarkersSigB_sorted.bam.bai"

#bowtie2 --very-sensitive-local -x SigBMarkers -U $S_aureus_readspath -S $outputpath$OutputSamFile

#samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

#samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

#samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

#OutputTSV="S_aureus_Coverage_SigBmarkers.tsv"

#python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

# S Lugdunensis 
#################
OutputSamFile="SLug_MarkersSigB.sam"
OutputBamFile="SLug_MarkersSigB.bam"
OutputBamFileSorted="SLug_MarkersSigB_sorted.bam"
OutputBamFileSortedBai="SLug_MarkersSigB_sorted.bam.bai"

bowtie2 --very-sensitive-local -x SigBMarkers -U $S_lug_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_Lugdunensis_Coverage_SigBmarkers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

# S capitis
#################
OutputSamFile="SCap_MarkersSigB.sam"
OutputBamFile="SCap_MarkersSigB.bam"
OutputBamFileSorted="SCap_MarkersSigB_sorted.bam"
OutputBamFileSortedBai="SCap_MarkersSigB_sorted.bam.bai"

bowtie2 --very-sensitive-local -x SigBMarkers -U $S_cap_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_capitis_Coverage_SigBmarkers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

# S haemolyticus
#################
OutputSamFile="SHaem_MarkersSigB.sam"
OutputBamFile="SHaem_MarkersSigB.bam"
OutputBamFileSorted="SHaem_MarkersSigB_sorted.bam"
OutputBamFileSortedBai="SHaem_MarkersSigB_sorted.bam.bai"

bowtie2 --very-sensitive-local -x SigBMarkers -U $S_haem_readspath -S $outputpath$OutputSamFile

samtools view -bS $outputpath$OutputSamFile > $outputpath$OutputBamFile

samtools sort  $outputpath$OutputBamFile >  $outputpath$OutputBamFileSorted

samtools index $outputpath$OutputBamFileSorted > $outputpath$OutputBamFileSortedBai

OutputTSV="S_haemolyticus_Coverage_SigBmarkers.tsv"

python /home/acampbe/pico_galaxy/tools/coverage_stats/coverage_stats.py -b $outputpath$OutputBamFileSorted -i $outputpath$OutputBamFileSortedBai -o $outputpath$OutputTSV

