# Amy Campbell
# March 2023
# Take metagenomic reads from Patient 176 and map them to genetic markers for CC1 and CC5 in Patient176

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv

# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/BTdb/"
mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient176/Patient176GeneMarkers_Final.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient176/output/"
mkdir -p $outputpath

# Build bowtie index of the markers
###################################
IndexName="Patient176GeneMarkers_Final"

bowtie2-build $markergenepath $IndexName
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment on 176_0 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_0.fastq.gz
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

# Run bowtie2 local alignment on 176_1 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_1.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_2 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_2.fastq.gz
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

# Run bowtie2 local alignment on 176_6 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_6.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_4 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_4.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_8 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_8.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_6 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_6.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_9 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_9.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_10 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_10.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_11 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_11.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_12 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_12.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

# Run bowtie2 local alignment on 176_13 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_176_13.fastq.gz
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
#samtools faidx $markergenepath

bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

BCFfolder="BCFs/"
mkdir -p $outputpath$BCFfolder
mv $outputpath*.bcf $outputpath$BCFfolder

# These intermediate files eat up so much storage
rm  $outputpath*.bam
rm  $outputpath*.sam

