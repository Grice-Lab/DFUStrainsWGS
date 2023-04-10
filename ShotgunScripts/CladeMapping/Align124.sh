# Amy Campbell
# March 2023
# Take metagenomic reads from Patient 124 and map them to genetic markers for CC1 and CC15 in Patient124

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv

# Remove any genes with >30x mean depth from 2million S. epi or S. pettenkefori reads
######################################################################################
oldmarkers="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/Patient124GeneMarkers.fasta"
SepicovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/Alignments/S_Epidermidis_Coverage_124markers.tsv"
SpetcovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/Alignments/S_Pettenkefori_Coverage_124markers.tsv"
Refgenome="DORN781"

Rscript RemoveTooConserved.R $oldmarkers $SepicovgStats $SpetcovgStats $Refgenome



# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/BTdb/"
mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/Patient124GeneMarkers_Final.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient124/output/"
mkdir -p $outputpath

# Build bowtie index of the markers
###################################
IndexName="Patient124GeneMarkers_Final"

bowtie2-build $markergenepath $IndexName
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment on 124_0 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_0.fastq.gz
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

rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext

# Run bowtie2 local alignment on 124_1 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_1.fastq.gz
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

#rm $outputpath$outputname$samext

rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext

# Run bowtie2 local alignment on 124_6 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_6.fastq.gz
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
rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext

#rm $outputpath$outputname$samext


# Run bowtie2 local alignment on 124_4 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_4.fastq.gz
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

#rm $outputpath$outputname$samext
rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext


# Run bowtie2 local alignment on 124_9 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_9.fastq.gz
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

#rm $outputpath$outputname$samext
rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext

# Run bowtie2 local alignment on 124_11 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_11.fastq.gz
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

#rm $outputpath$outputname$samext
rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext

# Run bowtie2 local alignment on 124_12 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_12.fastq.gz
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

#rm $outputpath$outputname$samext
rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext

# Run bowtie2 local alignment on 124_13 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_124_13.fastq.gz
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

#rm $outputpath$outputname$samext
rm $outputpath$outputname$samext
rm $outputpath$outputname$sortedbamext
rm $outputpath$outputname$bamext

