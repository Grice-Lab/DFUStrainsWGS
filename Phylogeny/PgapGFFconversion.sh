#!/bin/bash 
#Amy Campbell

filepath="/home/acampbe/DFU/data/WGS_2020/completed_PGAP"
outputpath="/home/acampbe/DFU/data/WGS_2020/ReformattedPGAP/"

mkdir -p $outputpath

for file in $filepath/*; do

	folder=$(basename $file)
	blank=""
	fextent="output"
	slash="/"
	genomename=${folder/$fextent/$blank}
	inputfilename="annot.gff"
	oldfile=$file$slash$inputfilename
	gffext=".gff"
	fastext="_pgap.fasta"

	newtemp=$outputpath$inputfilename
	newfile=$outputpath$genomename$gffext

	cp $oldfile $newtemp
	
	#sed -i 's/pgaptmp/'"$genomename"'/g' $newtemp	

	printf '%s\n' '##gff-version 3' > $newfile
	grep '##sequence-region' $newtemp >> $newfile
	grep "$(printf '\t')CDS" $newtemp >> $newfile
	printf '%s\n' '##FASTA' >> $newfile
	cat $file$slash$genomename$fastext >> $newfile
	rm $newtemp
done




#printf '%s\n' '##gff-version 3' > DORN244newtest.gff
#grep '##sequence-region' DORN244test.gff >> DORN244newtest.gff

#grep 'Protein Homology' DORN244test.gff >> DORN244newtest.gff
#printf '%s\n' '##FASTA' >> DORN244newtest.gff
#cat DORN244_pgap.fasta >> DORN244newtest.gff
