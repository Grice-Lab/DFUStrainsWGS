#!bin/bash

source ~/mambaforge/bin/activate DBSCANSWA

#contigs=/home/acampbe/DFU/data/WGS_2020/PhageResults/FixedFNAs/DORN925_fixed.fna
outputdir=/home/acampbe/DFU/data/WGS_2020/PhageResults/CheckVResults

mkdir -p $outputdir
export CHECKVDB=/home/acampbe/DownloadedDatabases/checkv-db-v1.5


for contigsfile in /home/acampbe/DFU/data/WGS_2020/PhageResults/FixedFNAs/* ; do
	basefilename=$(basename $contigsfile)
	fnaext="_fixed.fna"
	blank=""
	Prefix=${basefilename/$fnaext/$blank}
	outputfolder=$outputdir/$Prefix

	checkv end_to_end $contigsfile $outputfolder -t 16 
	
done


