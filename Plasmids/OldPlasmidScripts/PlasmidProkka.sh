# Amy Campbell
# 11/2020
# BLAST and Prokka on the circular plasmids extracted from each CC5 lineage

source /home/acampbe/software/miniconda3/bin/activate prokenv

mkdir -p /home/acampbe/CC5Plasmids/ProkkaResults


outdir="/home/acampbe/CC5Plasmids/ProkkaResults/"


for file in /home/acampbe/CC5Plasmids/*; do
	filebase=$(basename $file)	
        circext="_circular.fasta" 
        blank=""
        noext=${filebase/$circext/$blank}

	prokka --outdir $outdir$noext --force --prefix $noext --genus Staphylococcus $file
	
done


