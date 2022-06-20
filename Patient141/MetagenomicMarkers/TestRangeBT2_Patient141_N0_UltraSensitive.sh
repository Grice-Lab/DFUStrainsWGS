# Test combinations of the following parameters:

# Bowtie2 -i -i S,1,0.75, -i S,1,0.50 -i S,1,0.25, -i S,1,.01
# Bowtie2 -D :  (20, 25, 30)
# Bowtie2 -R :  (3,5,7,9)
# Bowtie2 -L : (20,18,16,14)
# BCFtools -Q : (13, 9, 5, 1)

# Ultra ultra sensitive
# D: 45
# R: 20
# L: 5
# Q: 1
# S2:.001


source /home/acampbe/software/miniconda3/bin/activate BowtieEnv


readspath="/home/acampbe/DFU/data/DFU_Metagenome_Microbes/Patient141/filtered_sorted_141-01.fastq"
readsname="filtered_sorted_141-01"

#/home/acampbe/DFU/data/DFU_Metagenome_Microbes/Patient141/*.fastq

underscore="_"

markerpath1="/home/acampbe/DFU/data/WGS_2020/Markers141/XanthinMarker1.fasta"
markerpath2="/home/acampbe/DFU/data/WGS_2020/Markers141/XanthinMarker2.fasta"
markerpath3="/home/acampbe/DFU/data/WGS_2020/Markers141/XanthinMarker3.fasta"

bowtiepath="/home/acampbe/DFU/data/WGS_2020/Markers141/BowtieDB/"
outputpath="/home/acampbe/DFU/data/WGS_2020/Markers141/BTalignments/N0UltraSensitive/"

mkdir -p $bowtiepath
mkdir -p $outputpath

export BOWTIE2_INDEXES=$bowtiepath

# Make bowtie DB for each marker
################################
bowtie2-build $markerpath1 XanthinMarker1
mv *.bt2 $bowtiepath

bowtie2-build $markerpath2 XanthinMarker2
mv *.bt2 $bowtiepath

bowtie2-build $markerpath3 XanthinMarker3
mv *.bt2 $bowtiepath

samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"
bcfext=".bcf"
samtools="samtools"


n0suffix="_N0"
ivalprefix="S,1,0."

for file in /home/acampbe/DFU/data/DFU_Metagenome_Microbes/Patient141/*.fastq; do

fastqext=".fastq"
blank=""
readspath=$file
readname=$(basename $file)
readsname=${readname/$fastqext/$blank}


# Marker 1
###########
dvalue=45
rvalue=20
lvalue=5
i2value="001"
qvalue=1
markername="XanthinMarker1"
ival=$ivalprefix$i2value    

iterationkey=$outputpath$readsname$underscore$markername$underscore$dvalue$underscore$rvalue$underscore$lvalue$underscore$i2value$underscore$qvalue$underscore$n0suffix
iterationkeyfinal=$iterationkey


bowtie2 -N 0 -D $dvalue -R $rvalue -L $lvalue -i $ival -x $markername -U $readspath -S $iterationkey$samext
samtools view -bS $iterationkey$samext > $iterationkey$bamext 

samtools sort $iterationkey$bamext > $iterationkey$sortedbamext
samtools mpileup -Q $qvalue -uf $filename $markerpath1 $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext

# Marker 2
###########
#dvalue=30
#rvalue=9
#lvalue=14
markername="XanthinMarker2"
#qvalue=1
#i2value="10"
#ival=$ivalprefix$i2value

iterationkey=$outputpath$readsname$underscore$markername$underscore$dvalue$underscore$rvalue$underscore$lvalue$underscore$i2value$underscore$qvalue$underscore$n0suffix
iterationkeyfinal=$iterationkey


bowtie2 -N 0 -D $dvalue -R $rvalue -L $lvalue -i $ival -x $markername -U $readspath -S $iterationkey$samext
samtools view -bS $iterationkey$samext > $iterationkey$bamext

samtools sort $iterationkey$bamext > $iterationkey$sortedbamext
samtools mpileup -Q $qvalue -uf $filename $markerpath2 $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext

# Marker 3
###########
#dvalue=30 
#rvalue=9
#lvalue=18
markername="XanthinMarker3"
#qvalue=5
#i2value="05"
ival=$ivalprefix$i2value

iterationkey=$outputpath$readsname$underscore$markername$underscore$dvalue$underscore$rvalue$underscore$lvalue$underscore$i2value$underscore$qvalue$underscore$n0suffix
iterationkeyfinal=$iterationkey

bowtie2 -N 0 -D $dvalue -R $rvalue -L $lvalue -i $ival -x $markername -U $readspath -S $iterationkey$samext
samtools view -bS $iterationkey$samext > $iterationkey$bamext

samtools sort $iterationkey$bamext > $iterationkey$sortedbamext
samtools mpileup -Q $qvalue -uf $filename $markerpath3 $iterationkey$sortedbamext | bcftools view -Ov -o $iterationkeyfinal$bcfext


done
