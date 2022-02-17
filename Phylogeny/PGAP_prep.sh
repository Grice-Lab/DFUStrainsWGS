
folderpath="/Users/amycampbell/Documents/PGAP/StaphGenomesAnnotation/"

for assembly in /Users/amycampbell/Documents/FinalSetDFUIsolates/*.fasta; do

	fname=$(basename $assembly)
	Oextension="_Final.fasta"
	noext=""
	genomename=${fname/$Oextension/$noext}
	slash="/"
	newext="_pgap.fasta"
	
	newfilename=$folderpath$genomename$slash

	mkdir -p $newfilename	
	cp $assembly $newfilename$genomename$newext

	

done 
