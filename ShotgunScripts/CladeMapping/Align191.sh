# Amy Campbell
# March 2023
# Take metagenomic reads from Patient 191 and map them to genetic markers for CC1 and CC15 in Patient191

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv

# Remove any genes with >30x mean depth from 2million S. epi or S. pettenkefori reads
######################################################################################
oldmarkers="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Patient191GeneMarkers.fasta"
SepicovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Alignments/S_Epidermidis_Coverage_191markers.tsv"
SpetcovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Alignments/S_Pettenkefori_Coverage_191markers.tsv"
Refgenome="DORN2205"

#Rscript RemoveTooConserved.R $oldmarkers $SepicovgStats $SpetcovgStats $Refgenome


# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/BTdb/"
mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/Patient191GeneMarkers_Final.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/output/"
mkdir -p $outputpath

# Build bowtie index of the markers
###################################
IndexName="Patient191GeneMarkers_Final"

bowtie2-build $markergenepath $IndexName
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment on 191_0 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_0.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

#bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

#samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
#samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
#samtools faidx $markergenepath

#bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 191_1 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_1.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

#bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

#samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
#samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
#samtools faidx $markergenepath

#bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext


# Run bowtie2 local alignment on 191_2 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_2.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

#bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

#samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
#samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
#samtools faidx $markergenepath

#bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 191_4 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_4.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

#bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

#samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
#samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
#samtools faidx $markergenepath

#bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext


# Run bowtie2 local alignment on 191_6 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_6.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 191_9 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_9.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 191_10 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_10.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext


# Run bowtie2 local alignment on 191_11 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_11.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext


# Run bowtie2 local alignment on 191_12 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_12.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 191_13 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_191_13.fastq.gz
ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
outputname=$metagenomename$markers

bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext





rm /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/output/*.sam
rm /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/output/*.bam
rm /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/output/*.bam.bai
