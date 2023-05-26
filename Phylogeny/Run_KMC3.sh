# amy campbell
# summarizing 31-length k-mers found in each 

source ~/mambaforgeOLD/bin/activate ~/mambaforgeOLD/envs/pankmer


fasta_ext=".fasta"
NewFolder="Kmers"
kmers_ext="_Kmers"
fastafolder="Fasta"
mkdir -p  /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/Kmers

for fasta in /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/Fasta/*; do 
	
	outputfile=${fasta/$fasta_ext/$kmers_ext}
	outputfile=${outputfile/$fastafolder/$NewFolder}

	kmc -k31 -fm $fasta $outputfile .

done







#kmc -k41 -fm /home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/Fasta/DORN286_Final.fasta 41mers .

