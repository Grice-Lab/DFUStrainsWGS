# Amy Campbell
# March 2023
# Take metagenomic reads from Patient 111 and map them to genetic markers for CC1 and CC15 in Patient111

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv

# Remove any genes with >30x mean depth from 2million S. epi or S. pettenkefori reads
######################################################################################
#oldmarkers="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/Patient111GeneMarkers.fasta"
#SepicovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/Alignments/S_Epidermidis_Coverage_111markers.tsv"
#SpetcovgStats="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/Alignments/S_Pettenkefori_Coverage_111markers.tsv"
#Refgenome="DORN152"

#Rscript RemoveTooConserved.R $oldmarkers $SepicovgStats $SpetcovgStats $Refgenome



# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/BTdb/"
#mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

markergenepath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/Patient111GeneMarkers_Final.fasta"
outputpath="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/output/"
#mkdir -p $outputpath

# Build bowtie index of the markers
###################################
IndexName="Patient111GeneMarkers_Final"

#bowtie2-build $markergenepath $IndexName
#mv *.bt2 $bowtiepath

# Run bowtie2 local alignment on 111_0 with very sensitive parameters
##############################################################################
#ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_111_0.fastq.gz
#ReadsFname=$(basename $ReadsPath)
#Readsext=".fastq.gz"
#blank=""
#metagenomename=${ReadsFname/$Readsext/$blank}

#samext=".sam"
#bamext=".bam"
#sortedbamext="_sorted.bam"
#bcfext=".bcf"
#markers="_marker"
#outputname=$metagenomename$markers

#bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext

#samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
#samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
#samtools faidx $markergenepath

#bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

#rm $outputpath$outputname$samext


# Run bowtie2 local alignment on 111_1 with very sensitive parameters
##############################################################################
#ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_111_1.fastq.gz
#ReadsFname=$(basename $ReadsPath)
#Readsext=".fastq.gz"
#blank=""
#metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
#bamext=".bam"
#sortedbamext="_sorted.bam"
#bcfext=".bcf"
#markers="_marker"
#outputname=$metagenomename$markers

#bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext
#samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
#samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
#samtools faidx $markergenepath

#bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

#rm $outputpath$outputname$samext

# Run bowtie2 local alignment on 111_2 with very sensitive parameters
##############################################################################
#ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_111_2.fastq.gz
#ReadsFname=$(basename $ReadsPath)
Readsext=".fastq.gz"
blank=""
#metagenomename=${ReadsFname/$Readsext/$blank}

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
markers="_marker"
#outputname=$metagenomename$markers

#bowtie2 -N 1 --very-sensitive-local -x $IndexName -U $ReadsPath -S $outputpath$outputname$samext
#samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
#samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
#samtools faidx $markergenepath

#bcftools mpileup --annotate FORMAT/DP4 -f $markergenepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext

#rm $outputpath$outputname$samext


# Run bowtie2 local alignment on 111_3 with very sensitive parameters
##############################################################################
ReadsPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/DFUMetagenomes/kalan01_DFUwgs_111_3.fastq.gz
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

