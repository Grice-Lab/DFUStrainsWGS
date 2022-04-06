# Amy Campbell
# March 2022
# Take simulated reads and map them to the markers 

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv


# Assign paths + make necessary folders
#######################################

bowtiepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB/"
markerspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/XanthinMarkers.fasta"
readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads"
#outputpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BTalignments/"
outputpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BTalignments/Markers/"

#/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads

mkdir -p /home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB
mkdir -p $outputpath

export BOWTIE2_INDEXES=$bowtiepath


# Build bowtie index
####################
bowtie2-build $markerspath XanthinMarkers
mv *.bt2 $bowtiepath


# Run bowtie2 local alignment on the simulated reads with default params for "very sensitive local"
####################################################################################################

for fwdreads in $readspath/*R1.fastq ; do 
	fwd=$fwdreads
	reverseext="_R2.fastq"
	fwdext="_R1.fastq"
	blank=""
	samext=".sam"
	bamext=".bam"
	sortedbamext="_sorted.bam"
	bcfext=".bcf"

	rev=${fwdreads/$fwdext/$reverseext}

	baseFilename=$(basename $fwdreads)
	genomename=${baseFilename/$fwdext/$blank}	

	#####################
	# Run local alignment
	#####################

	bowtie2 -N 1 --very-sensitive-local -x XanthinMarkers -1 $fwd -2 $rev -S $outputpath$genomename$samext

	# Convert to .bam file
	######################
	samtools view -bS $outputpath$genomename$samext > $outputpath$genomename$bamext	

	# Sort the bam file
	###################
	samtools sort $outputpath$genomename$bamext > $outputpath$genomename$sortedbamext

	# BCFTools to call variants
	###########################
	samtools faidx $markerspath
	bcftools mpileup -f $markerspath $outputpath$genomename$sortedbamext | bcftools view -Ov -o $outputpath$genomename$bcfext 

	# clean up intermediate files 
	rm $outputpath$genomename$samext
	rm $outputpath$genomename$sortedbamext
	rm $outputpath$genomename$bamext
done

