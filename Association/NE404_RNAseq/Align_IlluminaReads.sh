# Confirming Yabj/SpoVG are really absent from these genomes I think they're absent from 
source ~/mambaforge/bin/activate BowtieEnv23

bowtiepath="/home/acampbe/DFU/data/WGS_2020/yabJ_spoVG/BTdb/"
mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath


yabJ_spovG_marker="/home/acampbe/DFU/data/WGS_2020/yabJ_spoVG/yabJ_spoVGregion_USA300_FPR3757.fasta"
outputpath="/home/acampbe/DFU/data/WGS_2020/yabJ_spoVG/read_align_output/"
mkdir -p $outputpath
outputcounts="/home/acampbe/DFU/data/WGS_2020/yabJ_spoVG/read_align_output/NumBasesCovered10X.txt"

echo "start" > $outputcounts

# /project/grice/storage/HiSeq/WGS/HiSeq_19/



# Build bowtie index of the YabJSpoVG sequence from JE2
#######################################################
IndexName="yabJspoVG"

bowtie2-build $yabJ_spovG_marker $IndexName
mv *.bt2 $bowtiepath


# DORN1499 (negative)
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1499trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1499trimmedgalore_val_2.fastq

# Run bowtie2
#############
basenamefile=$(basename $readspath1)
trimmed1ext="trimmedgalore_val_1.fastq"
blank=""
noext=${basenamefile/$trimmed1ext/$blank}

samext=".sam"
bamext=".bam"


bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment
#####################################################################
totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext >> $outputcounts
echo $totalBasesCovered10x >> $outputcounts


# DORN1743 (negative)
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1743trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1743trimmedgalore_val_2.fastq

# Run bowtie2
#############
basenamefile=$(basename $readspath1)
trimmed1ext="trimmedgalore_val_1.fastq"
blank=""
noext=${basenamefile/$trimmed1ext/$blank}

samext=".sam"
bamext=".bam"


bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment
#####################################################################
totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext >> $outputcounts
echo $totalBasesCovered10x >> $outputcounts


# DORN1761 (negative)
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1761trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1761trimmedgalore_val_2.fastq

# Run bowtie2
#############
basenamefile=$(basename $readspath1)
trimmed1ext="trimmedgalore_val_1.fastq"
blank=""
noext=${basenamefile/$trimmed1ext/$blank}

samext=".sam"
bamext=".bam"


bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment
#####################################################################
totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext >> $outputcounts
echo $totalBasesCovered10x >> $outputcounts


# DORN1952 (negative)
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1952trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1952trimmedgalore_val_2.fastq

# Run bowtie2
#############
basenamefile=$(basename $readspath1)
trimmed1ext="trimmedgalore_val_1.fastq"
blank=""
noext=${basenamefile/$trimmed1ext/$blank}

samext=".sam"
bamext=".bam"


bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment
#####################################################################
totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext >> $outputcounts
echo $totalBasesCovered10x >> $outputcounts


# DORN2221 (negative)
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2221trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2221trimmedgalore_val_2.fastq

# Run bowtie2
#############
basenamefile=$(basename $readspath1)
trimmed1ext="trimmedgalore_val_1.fastq"
blank=""
noext=${basenamefile/$trimmed1ext/$blank}

samext=".sam"
bamext=".bam"


bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment
#####################################################################
totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext >> $outputcounts
echo $totalBasesCovered10x >> $outputcounts


# DORN1729 (positive CC8 case)
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1729trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1729trimmedgalore_val_2.fastq

# Run bowtie2
#############
basenamefile=$(basename $readspath1)
trimmed1ext="trimmedgalore_val_1.fastq"
blank=""
noext=${basenamefile/$trimmed1ext/$blank}

samext=".sam"
bamext=".bam"


bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment
#####################################################################
totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext >> $outputcounts
echo $totalBasesCovered10x >> $outputcounts


# DORN1869 (positive CC5 case)
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1869trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1869trimmedgalore_val_2.fastq

# Run bowtie2
#############
basenamefile=$(basename $readspath1)
trimmed1ext="trimmedgalore_val_1.fastq"
blank=""
noext=${basenamefile/$trimmed1ext/$blank}

samext=".sam"
bamext=".bam"


bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment
#####################################################################
totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
echo $noext >> $outputcounts
echo $totalBasesCovered10x >> $outputcounts


# Run bowtie2
#############
#basenamefile=$(basename $readspath1)
#trimmed1ext="trimmedgalore_val_1.fastq"
#blank=""
#noext=${basenamefile/$trimmed1ext/$blank}

#samext=".sam"
#bamext=".bam"


#bowtie2 -x $IndexName -1 $readspath1 -2 $readspath2 -S $outputpath$noext$samext

#samtools sort $outputpath$noext$samext > $outputpath$noext$bamext
#samtools index $outputpath$noext$bamext

# Get number of bases in reference with >=10x coverage by the alignment 
#####################################################################
#totalBasesCovered10x=$(samtools mpileup $outputpath$noext$bamext | awk -v X=9 '$4>x' | wc -l)
#echo $noext >> $outputcounts
#echo $totalBasesCovered10x >> $outputcounts


