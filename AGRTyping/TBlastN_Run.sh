#!/bin/bash
# Blasting each AGR type against the contigs (uncleaned) 
source /home/acampbe/software/miniconda3/bin/activate BlastEnv


for fastaname in /home/acampbe/DFU/data/WGS_2020/BlastP_AGR_D/databases/*.fasta; do
	foldername=$(basename $fastaname)
	ext="_contigs.fasta"
	replace=""
	nameToUse=${foldername/$ext/$replace}
	marker=">"
	delimiter="_"
	sed -i "s/$marker/$marker$nameToUse$delimiter/g" $fastaname
done

for fastaname in /home/acampbe/DFU/data/WGS_2020/BlastP_AGR_D/databases/*.fasta; do
	foldername=$(basename $fastaname)
	makeblastdb -in $fastaname -input_type fasta -dbtype nucl -out $foldername 

done
 
#for fastaname in /home/acampbe/DFU/data/WGS_2020/BlastP_AGR_D/controls/*.fasta; do
#        foldername=$(basename $fastaname)
#        makeblastdb -in $fastaname -input_type fasta -dbtype nucl -out $foldername

#done
for fastaname in /home/acampbe/DFU/data/WGS_2020/BlastP_AGR_D/databases/*.fasta; do
	base=$(basename $fastaname)
	ext=".txt"
	tblastn -query /home/acampbe/DFU/data/WGS_2020/BlastP_AGR_D/AGR_D_aa.faa -out $base$ext -db $base -outfmt "6 qseqid sseqid pident qcovs evalue bitscore"

	
done

for fastaname in /home/acampbe/DFU/data/WGS_2020/BlastP_AGR_D/controls/*.fasta; do
        base=$(basename $fastaname)
        ext=".txt"
        tblastn -query /home/acampbe/DFU/data/WGS_2020/BlastP_AGR_D/AGR_D_aa.faa -out $base$ext -db $base -outfmt "6 qseqid sseqid pident qcovs evalue bitscore"


done
