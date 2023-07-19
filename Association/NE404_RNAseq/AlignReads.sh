# Amy Campbell
# Align RNAseq paired end, trimmed reads of JE2 vs. NE404
#  to reference for FPR3757
# 07/23

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv


bowtiepath="/home/acampbe/DFU/data/RNASeqNe404/BTindices"
readspath="/project/grice/storage/RNASeq_JE2_NE404/Trimmed/"

outputpath="/home/acampbe/DFU/data/RNASeqNe404/BTAlign/"
refpath="/home/acampbe/DFU/data/RNASeqNe404/USA300_FPR3757.fasta"

mkdir -p $bowtiepath
mkdir -p $outputpath

export BOWTIE2_INDEXES=$bowtiepath

#/project/grice/storage/RNASeq_JE2_NE404/Trimmed/JE2_rep3_R1_trimmed.fq.gz

# Build bowtie index
####################
refname="FPR3757ref"

bowtie2-build $refpath $refname
mv *.bt2 $bowtiepath

forwardext="_val_1.fq.gz"
revext="_val_2.fq.gz"
noext=""
samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"

for filename in $readspath/*_val_1.fq.gz; do 
	base=$(basename $filename)
	prefixname=${base/$forwardext/$noext}
	
	forwardfile=$filename
	reversefile=${filename/$forwardext/$revext}

	bowtie2 -x $refname -1 $forwardfile -2 $reversefile -S $outputpath$prefixname$samext
	samtools view -bS $outputpath$prefixname$samext > $outputpath$prefixname$bamext	
	samtools sort $outputpath$prefixname$bamext > $outputpath$prefixname$sortedbamext
	rm $outputpath$prefixname$samext

done

