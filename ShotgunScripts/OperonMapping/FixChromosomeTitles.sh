# To make the publicly-downloaded USA300_FPR3757 genome compatible with formatting
# in the SNPeff prebuilt database, change the following chromosomes' names:
# NC_007793.1 --> Chromosome
# NC_007790.1 --> pUSA01 
# NC_007791.1 --> pUSA02
# NC_007792.1 --> pUSA03

for file in /home/acampbe/DFU/data/StressOperonMarkers/Alignments/FullGenome_output/*.vcf ; do 
	sed -i 's/NC_007793.1/Chromosome/g' $file
	sed -i 's/NC_007790.1/pUSA01/g' $file
	sed -i 's/NC_007791.1/pUSA02/g' $file
	sed -i 's/NC_007792.1/pUSA03/g' $file
done
