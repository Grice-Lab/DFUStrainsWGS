#!bin/bash
# Amy Campbell April 2023
# Take in a metagenomic sample and align it to full USA300/LAC strain NCBI genome
# To make a VCF (deleting intermediates) for input into SnpEff

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomes/StressOperonMarkers/BTdb/"
export BOWTIE2_INDEXES=$bowtiepath

markergenepath="/home/acampbe/software/snpEff/data/CP000255.1/USA300_LAC.fasta"
outputpath="/home/acampbe/DFU/data/StressOperonMarkers/Alignments/FullGenome_output/"

#inputreadspath=$1
#file=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_171_2.fastq.gz


for file in /home/acampbe/DFU/data/DFU_Metagenome_Microbes/*.fastq.gz; do

	Readspath=$file
        OutputPathBase=$(basename $Readspath)
        fastq=".fastq.gz"
        blank=""
        OutputPrefix=${OutputPathBase/$fastq/$blank}
        OutputFilePrefix=$outputpath$OutputPrefix
        
        samext="_LAC.sam"
        bamext="_LAC.bam"
        sortedbamext="_LAC_sorted.bam"
        bcfext="_LAC.vcf"

        bowtie2 -N 1 --very-sensitive-local -x LACref -U $Readspath -S $OutputFilePrefix$samext
        samtools view -bS $OutputFilePrefix$samext > $OutputFilePrefix$bamext
        samtools sort $OutputFilePrefix$bamext > $OutputFilePrefix$sortedbamext
	samtools faidx $markergenepath
	bcftools mpileup -f $markergenepath $OutputFilePrefix$sortedbamext | bcftools view -Ov -o $OutputFilePrefix$bcfext

#        bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $OutputFilePrefix$sortedbamext | bcftools view -Ov -o $OutputFilePrefix$bcfext

        rm $OutputFilePrefix$samext
        rm $OutputFilePrefix$bamext
        rm $OutputFilePrefix$sortedbamext
	echo "Complete"
	echo $OutputPathBase
done
