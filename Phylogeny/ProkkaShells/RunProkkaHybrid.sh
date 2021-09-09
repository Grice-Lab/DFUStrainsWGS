#!/bin/bash
mkdir -p /home/acampbe/DFU/data/WGS_2020/ONP/Prokka
source /home/acampbe/software/miniconda3/bin/activate prokenv

fastaext=".fasta"
blank=""
outputdir="/home/acampbe/DFU/data/WGS_2020/ONP/Prokka/"

for fast in /home/acampbe/DFU/data/WGS_2020/HybridAssemblies/*; do
	assemblyname=$(basename $fast)
	genomename=${assemblyname/$fastaext/$blank}
	prokka --outdir $outputdir$genomename --force --prefix $genomename --genus Staphylococcus $fast
	echo $genomename
done


