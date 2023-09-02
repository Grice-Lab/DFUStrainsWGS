#!/bin/bash
# Get mash distances between every identified phage
# so that we can look at the mash distance between members of the same phage cluster

source ~/mambaforge/bin/activate DBSCANSWA


allfastas=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/IndividualFastas
outputfileMashSketch=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/Phages.msh
outputfolderTSVs=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/mashTSVs/
outputTSV="/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults_Parsed/cdhit_80/MashDistsAll1000.tsv"

echo "phageinstance1	phageinstance2	distance	p	NumberHashes" > $outputTSV

mkdir -p $outputfolderTSVs

# Make a mash sketch of all the phages found
#############################################
mash sketch -k 21 -s 1000 $allfastas/*.fasta -o $outputfileMashSketch

fastaext=".fasta"
newext="_dists1000.tsv"
for fastaitem in $allfastas/*.fasta ; do 
	basestring=$(basename $fastaitem)
	tsvfilename=${basestring/$fastaext/$newext}
	mash dist $outputfileMashSketch $fastaitem > $outputfolderTSVs$tsvfilename
	cat $outputfolderTSVs$tsvfilename >> $outputTSV
	rm $outputfolderTSVs$tsvfilename


done
