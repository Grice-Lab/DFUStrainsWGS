# Amy Campbell
#  Step 2 in the plan laid out in 
# "2022-04-03 Testing range of sensitivity
# parameters for staphyloxanthin markers"
# Using very-sensitive-local Bowtie2 parameters,
# align the sorted/combined simulated shotgun reads
# for each genome

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv

# DEFINE PATHS
##############

reads925="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/DORN925sorted.fastq"
reads1088="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/DORN1088sorted.fastq"

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


# BUILD BOWTIE INDICES
######################

bowtie2-build $genomepath $genomename
bowtie2-build $markerpath $markername

mv *.bt2 $bowtiepath


# ALIGN EACH TO REFERENCE 925
#############################

# DORN925
#readsname925="DORN925Reads_"

#bowtie2 -N 1 --very-sensitive-local -x $genomename -U $reads925 -S $outputpath$readsname925$genomename$samext

#samtools view -bS $outputpath$readsname925$genomename$samext > $outputpath$readsname925$genomename$bamext

#samtools sort $outputpath$readsname925$genomename$bamext > $outputpath$readsname925$genomename$sortedbamext
      
#samtools faidx $genomepath

#bcftools mpileup -f $genomepath $outputpath$readsname925$genomename$sortedbamext | bcftools view -Ov -o $outputpath$readsname925$genomename$bcfext


# DORN1088 reads

readsname1088="DORN1088Reads_"

bowtie2 -N 1 --very-sensitive-local -x $genomename -U $reads1088 -S $outputpath$readsname1088$genomename$samext

samtools view -bS $outputpath$readsname1088$genomename$samext > $outputpath$readsname1088$genomename$bamext

samtools sort $outputpath$readsname1088$genomename$bamext > $outputpath$readsname1088$genomename$sortedbamext

bcftools mpileup -f $genomepath $outputpath$readsname1088$genomename$sortedbamext | bcftools view -Ov -o $outputpath$readsname1088$genomename$bcfext


# ALIGN EACH TO THE 500 BASE MARKERS
####################################

# DORN1088
readsname1088="DORN1088Reads_"

bowtie2 -N 1 --very-sensitive-local -x $markername -U $reads1088 -S $outputpath$readsname1088$markername$samext

samtools view -bS $outputpath$readsname1088$markername$samext > $outputpath$readsname1088$markername$bamext

samtools sort $outputpath$readsname1088$markername$bamext > $outputpath$readsname1088$markername$sortedbamext

samtools faidx $markerpath

bcftools mpileup -f $markerpath $outputpath$readsname1088$markername$sortedbamext | bcftools view -Ov -o $outputpath$readsname1088$markername$bcfext





# DORN925
readsname925="DORN925Reads_"

bowtie2 -N 1 --very-sensitive-local -x $markername -U $reads925 -S $outputpath$readsname925$markername$samext

samtools view -bS $outputpath$readsname925$markername$samext > $outputpath$readsname925$markername$bamext

samtools sort $outputpath$readsname925$markername$bamext > $outputpath$readsname925$markername$sortedbamext

samtools faidx $markerpath

bcftools mpileup -f $markerpath $outputpath$readsname925$markername$sortedbamext | bcftools view -Ov -o $outputpath$readsname925$markername$bcfext



rm $outputpath$readsname1088$markername$sortedbamext
rm $outputpath$readsname1088$markername$samext
rm $outputpath$readsname1088$markername$bamext


rm $outputpath$readsname925$markername$sortedbamext
rm $outputpath$readsname925$markername$samext
rm $outputpath$readsname925$markername$bamext


rm $outputpath$readsname1088$genomename$sortedbamext
rm $outputpath$readsname1088$genomename$samext
rm $outputpath$readsname1088$genomename$bamext


rm $outputpath$readsname925$genomename$sortedbamext
rm $outputpath$readsname925$genomename$samext
rm $outputpath$readsname925$genomename$bamext


