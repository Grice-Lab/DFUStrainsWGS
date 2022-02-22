for folder in /home/acampbe/DFU/data/WGS_2020/completed_prokka/gff_files/* ; do 
	slash="/"
	fnamext=".gff"
	base=$(basename $folder)	
	cp $folder$slash$base$fnamext /home/acampbe/DFU/data/WGS_2020/FinalProkka/WithEpi/
	cp $folder$slash$base$fnamext /home/acampbe/DFU/data/WGS_2020/FinalProkka/WithoutEpi/

done

rm /home/acampbe/DFU/data/WGS_2020/FinalProkka/WithoutEpi/S_epidermidis.gff
