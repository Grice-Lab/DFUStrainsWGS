#!/bin/bash
# Run DBSCAN-SWA on each of the 220 genomes

source ~/mambaforge/bin/activate DBSCANSWA

fastasuffix="_Final.fasta"
blank=""
outputfolderPrefix="/home/acampbe/DFU/data/WGS_2020/PhageResults/"


for fastafile in /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/*.fasta ; do

	baseFilename=$(basename $fastafile)
	otptprefix=${baseFilename/$fastasuffix/$blank}
	#echo $fastafile
	#echo $outputfolderPrefix$otptprefix
	#echo $otptprefix
	
	python3 /home/acampbe/DBSCAN-SWA/bin/dbscan-swa.py --input $fastafile --output $outputfolderPrefix$otptprefix --prefix $otptprefix

done
