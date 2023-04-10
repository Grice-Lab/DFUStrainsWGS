# Amy Campbell
# March 2023
# Take metagenomic reads from Patient 159 and map them to genetic markers for CC1 and CC15 in Patient159

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv

# Remove any genes with >30x mean depth from 2million S. epi or S. pettenkefori reads
######################################################################################
oldmarkers="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159/Patient159GeneMarkers.fasta"
SepicovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159/Alignments/S_Epidermidis_Coverage_159markers.tsv"
SpetcovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159/Alignments/S_Pettenkefori_Coverage_159markers.tsv"
Refgenome="DORN1364"

Rscript RemoveTooConserved.R $oldmarkers $SepicovgStats $SpetcovgStats $Refgenome



# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/BTdb/"
mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159/Patient159GeneMarkers_Final.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159/output/"
mkdir -p $outputpath

# Build bowtie index of the markers
###################################
IndexName="Patient159GeneMarkers_Final"

bowtie2-build $markergenepath $IndexName
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment on 159_0 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_159_0.fastq.gz
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

# Run bowtie2 local alignment on 159_1 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_159_1.fastq.gz
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


# Run bowtie2 local alignment on 159_2 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_159_2.fastq.gz
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

# Run bowtie2 local alignment on 159_4 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_159_4.fastq.gz
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

