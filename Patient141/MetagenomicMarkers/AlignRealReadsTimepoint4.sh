# Amy Campbell
# March 2022
# Take simulated reads and map them to the markers 

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv

# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB/"
markerspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/XanthinMarkers.fasta"
genomepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/DORN925_Final.fasta"
#readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads"
readspath="/home/acampbe/DFU/data/DFU_Metagenome_Microbes/Patient141/filtered_sorted_141-04.fastq"
outputpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/MetagenomeRealAlignments/"

#/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads

mkdir -p /home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB
mkdir -p $outputpath

export BOWTIE2_INDEXES=$bowtiepath

#/home/acampbe/DFU/data/DFU_Metagenome_Microbes/Timepoint01/filtered_sorted_141-10.fastq.gz


# Build bowtie index
####################
genomename="DORN925_Final"

bowtie2-build $genomepath $genomename
mv *.bt2 $bowtiepath

# Run bowtie2 local alignment on the real reads with very sensitive parameters
##############################################################################

metagenomename="filtered_sorted_141-04"
samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
genomename="DORN925_Final"
outputname=$genomename$metagenomename
markername="XanthinMarkers"
outputname2=$markername$metagenomename

bowtie2 -N 1 --very-sensitive-local -x $genomename -U $readspath -S $outputpath$outputname$samext

samtools view -bS $outputpath$outputname$samext > $outputpath$outputname$bamext
samtools sort $outputpath$outputname$bamext > $outputpath$outputname$sortedbamext
samtools faidx $genomename
bcftools mpileup -f $genomepath $outputpath$outputname$sortedbamext | bcftools view -Ov -o $outputpath$outputname$bcfext


# Just markers
bowtie2-build $markerspath $markername
mv *.bt2 $bowtiepath
bowtie2 -N 1 --very-sensitive-local -x $markername -U $readspath -S $outputpath$outputname2$samext
samtools view -bS $outputpath$outputname2$samext > $outputpath$outputname2$bamext
samtools sort $outputpath$outputname2$bamext > $outputpath$outputname2$sortedbamext
samtools faidx $markername
bcftools mpileup -f $markerspath $outputpath$outputname2$sortedbamext | bcftools view -Ov -o $outputpath$outputname2$bcfext





# Run bowtie2 local alignment on the simulated reads with default params
########################################################################

#for fwdreads in $readspath/*R1.fastq ; do 
#	fwd=$fwdreads
#	reverseext="_R2.fastq"
#	fwdext="_R1.fastq"
#	blank=""
#	samext=".sam"
#	bamext=".bam"
#	sortedbamext="_sorted.bam"
#	bcfext=".bcf"

#	rev=${fwdreads/$fwdext/$reverseext}

#	baseFilename=$(basename $fwdreads)
#	genomename=${baseFilename/$fwdext/$blank}	

	#####################
	# Run local alignment
	#####################

#	bowtie2 -N 1 --very-sensitive-local -x XanthinMarkers -1 $fwd -2 $rev -S $outputpath$genomename$samext

	# Convert to .bam file
	######################
#	samtools view -bS $outputpath$genomename$samext > $outputpath$genomename$bamext	

	# Sort the bam file
	###################
#	samtools sort $outputpath$genomename$bamext > $outputpath$genomename$sortedbamext

	# BCFTools to call variants
	###########################
#	samtools faidx $markerspath
#	bcftools mpileup -f $markerspath $outputpath$genomename$sortedbamext | bcftools view -Ov -o $outputpath$genomename$bcfext 

#done

