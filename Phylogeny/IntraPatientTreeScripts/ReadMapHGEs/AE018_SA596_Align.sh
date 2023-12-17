source ~/mambaforge/bin/activate BowtieEnv23

bowtiepath=/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/BTdbs/

mkdir -p $bowtiepath

export BOWTIE2_INDEXES=$bowtiepath



outputcounts=AE018_SA596_Coverage.txt
touch outputcounts
trimmed1ext="trimmedgalore_val_1.fastq"

samext=".sam"

bamext=".bam"

blank=""

IndexName="AE018_SA596"

markerpath=/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/Markers/AE018_SA596.fasta
bowtie2-build $markerpath $IndexName
mv *.bt2 $bowtiepath


readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN596trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN596trimmedgalore_val_1.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN620trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN620trimmedgalore_val_1.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN613trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN613trimmedgalore_val_1.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


