readsfile1="/home/acampbe/DFU/data/DFU_Metagenome_Microbes/Patient141/filtered_sorted_141-01.fastq"
readsname1="filtered_sorted_141-01"

readsfile2="/home/acampbe/DFU/data/DFU_Metagenome_Microbes/Patient141/filtered_sorted_141-04.fastq"
readsname2="filtered_sorted_141-04"

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv


genomepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/DORN925_Final.fasta"
markerpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/XanthinMarkersExtended.fasta"
genomename="DORN925Full"
markername="XanthinMarkersExt"

bowtiepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB/"

outputpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BTalignments/TestMarkerVsGenome/"

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"

mkdir -p $outputpath

export BOWTIE2_INDEXES=$bowtiepath



bowtie2 -N 1 --very-sensitive-local -x $markername -U $readsfile2 -S $outputpath$readsname2$markername$samext

samtools view -bS $outputpath$readsname2$markername$samext > $outputpath$readsname2$markername$bamext

samtools sort $outputpath$readsname2$markername$bamext > $outputpath$readsname2$markername$sortedbamext

samtools faidx $markerpath

bcftools mpileup -f $markerpath $outputpath$readsname2$markername$sortedbamext | bcftools view -Ov -o $outputpath$readsname2$markername$bcfext

