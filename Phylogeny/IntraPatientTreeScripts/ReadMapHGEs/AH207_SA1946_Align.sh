source ~/mambaforge/bin/activate BowtieEnv23

bowtiepath=/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/BTdbs/

mkdir -p $bowtiepath

export BOWTIE2_INDEXES=$bowtiepath



outputcounts=AH207_SA1946_Coverage.txt
touch $soutputcounts
trimmed1ext="trimmedgalore_val_1.fastq"

samext=".sam"

bamext=".bam"

blank=""

IndexName="AH207_SA1946"

markerpath=/home/acampbe/DFU/data/WGS_2020/PhagePlasmidMapping/Markers/AH207_SA1946.fasta
bowtie2-build $markerpath $IndexName
mv *.bt2 $bowtiepath


readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1946trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1946trimmedgalore_val_2.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


rm $outputpath$noext$samext
rm $outputpath$noext$bamext
rm $outputpath$noext$bamext.bai
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2221trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2221trimmedgalore_val_2.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


rm $outputpath$noext$samext
rm $outputpath$noext$bamext
rm $outputpath$noext$bamext.bai
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2130trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2130trimmedgalore_val_2.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


rm $outputpath$noext$samext
rm $outputpath$noext$bamext
rm $outputpath$noext$bamext.bai
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2093trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2093trimmedgalore_val_2.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


rm $outputpath$noext$samext
rm $outputpath$noext$bamext
rm $outputpath$noext$bamext.bai
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2058trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2058trimmedgalore_val_2.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


rm $outputpath$noext$samext
rm $outputpath$noext$bamext
rm $outputpath$noext$bamext.bai
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1961trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1961trimmedgalore_val_2.fastq
basenamefile=$(basename $readspath1)
noext=${basenamefile/$trimmed1ext/$blank}
bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext
samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext	$totalBasesCovered10x >> $outputcounts


rm $outputpath$noext$samext
rm $outputpath$noext$bamext
rm $outputpath$noext$bamext.bai
