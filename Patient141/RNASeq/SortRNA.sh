source /home/acampbe/miniconda3/bin/activate Patient141New

readspath1="/home/acampbe/DFU/data/RNASeq/Trimmed/ARM1trimmedBACKUP_R1.fastq"
readspath2="/home/acampbe/DFU/data/RNASeq/Trimmed/ARM1trimmedBACKUP_R2.fastq"

otputfwd="/home/acampbe/DFU/data/RNASeq/Trimmed/ARM1sortedBACKUP_R1"
otputrev="/home/acampbe/DFU/data/RNASeq/Trimmed/ARM1sortedBACKUP_R2"

#refpath1="/home/acampbe/DownloadedDatabases/SortMeDBs/rfam-5.8s-database-id98.fasta,rfam-5.8s-database-id98-db"
#refpath2="/home/acampbe/DownloadedDatabases/SortMeDBs/rfam-5s-database-id98.fasta,rfam-5s-database-id98-db"
#refpath3="/home/acampbe/DownloadedDatabases/SortMeDBs/silva-bac-16s-id90.fasta,silva-bac-16s-id90-db"
#refpath4="/home/acampbe/DownloadedDatabases/SortMeDBs/silva-bac-23s-id98.fasta,silva-bac-23s-id98-db"

refpath1="/home/acampbe/DownloadedDatabases/SortMeDBs/rfam-5.8s-database-id98.fasta"
refpath2="/home/acampbe/DownloadedDatabases/SortMeDBs/rfam-5s-database-id98.fasta"
refpath3="/home/acampbe/DownloadedDatabases/SortMeDBs/silva-bac-16s-id90.fasta"
refpath4="/home/acampbe/DownloadedDatabases/SortMeDBs/silva-bac-23s-id98.fasta"


/home/acampbe/software/sortmerna-4.3.4-Linux/bin/sortmerna --fastx --ref $refpath1 --ref $refpath2 --ref $refpath3 --ref $refpath4 --reads $readspath1  --reads $readspath2
#/home/acampbe/software/sortmerna-4.3.4-Linux/bin/sortmerna
