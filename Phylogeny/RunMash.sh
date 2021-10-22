source /home/acampbe/software/miniconda3/bin/activate TreeEnv

outputfolder="/home/acampbe/Mash/MashOutput/"
for f in /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Assemblies/*.fasta ; do
	inputname=$(basename $f)
	output="_Mash.tab"
	inputextension=".fasta"
	
	outputfilename=${inputname/$inputextension/$output}
	outputpath=$outputfolder$outputfilename
	
	mash screen -w -p 4 refseq.genomes.k21s1000.msh $f  > $outputpath

done
