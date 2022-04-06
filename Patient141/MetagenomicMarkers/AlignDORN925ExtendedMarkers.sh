# Testing that the reads are what I think they are 
source /home/acampbe/software/miniconda3/bin/activate BowtieEnv

# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB/"
markerspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/XanthinMarkers.fasta"
#readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads"
genomepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/XanthinMarkersExtended.fasta"

fwdpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/DORN925_R1.fastq"
revpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/DORN925_R2.fastq"

#outputpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/MetagenomeRealAlignments/"
outputpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BTalignments/Markers/"

mkdir -p $outputpath

#/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads

mkdir -p /home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB
#mkdir -p $outputpath

export BOWTIE2_INDEXES=$bowtiepath

#/home/acampbe/DFU/data/DFU_Metagenome_Microbes/Timepoint01/filtered_sorted_141-10.fastq.gz


# Build bowtie index
####################
genomename="DORN925XanthinMarkersExtended"

bowtie2-build $genomepath $genomename
mv *.bt2 $bowtiepath



# Run bowtie2 local alignment on the real reads with very sensitive parameters
##############################################################################

genomename="DORN925XanthinMarkersExtended"
samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"

#MIN_COVERAGE_DEPTH5=5
#MIN_COVERAGE_DEPTH20=20
#MIN_COVERAGE_DEPTH50=50

bowtie2 -N 1 --very-sensitive-local -x $genomename -1 $fwdpath -2 $revpath -S $outputpath$genomename$samext

samtools view -bS $outputpath$genomename$samext > $outputpath$genomename$bamext

samtools sort $outputpath$genomename$bamext > $outputpath$genomename$sortedbamext
      
samtools faidx $genomepath

bcftools mpileup -f $genomepath $outputpath$genomename$sortedbamext | bcftools view -Ov -o $outputpath$genomename$bcfext

# clean up intermediate files
rm $outputpath$genomename$samext
rm $outputpath$genomename$sortedbamext
rm $outputpath$genomename$bamext
