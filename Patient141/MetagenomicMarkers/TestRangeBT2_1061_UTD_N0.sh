# Test combinations of the following parameters:

# Bowtie2 -i -i S,1,0.75, -i S,1,0.50 -i S,1,0.25, -i S,1,.01
# Bowtie2 -D :  (20, 25, 30)
# Bowtie2 -R :  (3,5,7,9)
# Bowtie2 -L : (20,18,16,14)
# BCFtools -Q : (13, 9, 5, 1)

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv

readspath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads/DORN1061sorted.fastq"
readsname="DORN1061reads"

underscore="_"

markerpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/XanthinMarkersExtended.fasta"
markername="XanthinMarkersExt"
bowtiepath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BowtieDB/"
outputpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/BTalignments/TestRangeParameters/"

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
samtools="samtools"
mkdir -p $outputpath
samtoolsd="samtoolsdefault"

export BOWTIE2_INDEXES=$bowtiepath

iterationkey1=$outputpath$readsname$underscore$samtoolsd

#bowtie2 -N 1 --very-sensitive-local -x $markername -U $readspath -S $iterationkey1$samext
#samtools view -bS $iterationkey1$samext > $iterationkey1$bamext
#samtools sort $iterationkey1$bamext > $iterationkey1$sortedbamext
#samtools mpileup -uf $filename $markerpath $iterationkey1$sortedbamext | bcftools view -Ov -o $iterationkey1$bcfext

ivalprefix="S,1,0."
n0suffix="_N0"
for dvalue in 20 25 30; do

  for rvalue in 3 5 7 9; do

    for lvalue in 20 18 16 14; do

	for i2value in "75" "50" "25" "10" "05" ; do 

		ival=$ivalprefix$i2value    
        	iterationkey=$outputpath$readsname$underscore$dvalue$underscore$rvalue$underscore$lvalue$underscore$i2value$n0suffix


        	bowtie2 -N 1 -D $dvalue -R $rvalue -L $lvalue -i $ival -x $markername -U $readspath -S $iterationkey$samext
        	samtools view -bS $iterationkey$samext > $iterationkey$bamext 
        	samtools sort $iterationkey$bamext > $iterationkey$sortedbamext


		qvalue=13
        	iterationkeyfinal=$iterationkey$underscore$qvalue
        	samtools mpileup -Q $qvalue -uf $filename $markerpath $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext
	
        	qvalue=9
        	iterationkeyfinal=$iterationkey$underscore$qvalue
        	samtools mpileup -Q $qvalue -uf $filename $markerpath $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext

        	qvalue=5
        	iterationkeyfinal=$iterationkey$underscore$qvalue
        	samtools mpileup -Q $qvalue -uf $filename $markerpath $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext

		qvalue=1
        	iterationkeyfinal=$iterationkey$underscore$qvalue
        	samtools mpileup -Q $qvalue -uf $filename $markerpath $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext

		qvalue=0
		iterationkeyfinal=$iterationkey$underscore$qvalue
        	samtools mpileup -Q $qvalue -uf $filename $markerpath $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext


	
		rm $iterationkey$samext
		rm $iterationkey$bamext
		rm $iterationkey$sortedbamext




	done
    done

  done

done



