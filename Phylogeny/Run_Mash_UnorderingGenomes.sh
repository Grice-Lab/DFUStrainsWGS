# Screening two DORN isolates (presumed S. aureus) whose contigs failed to reorder against a S. aureus reference genome

source /home/acampbe/software/miniconda3/bin/activate TreeEnv

mash screen -w -p 4 refseq.genomes.k21s1000.msh /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Assemblies/DORN685_contigs.fasta
mash screen -w -p 4 refseq.genomes.k21s1000.msh /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Assemblies/DORN946_contigs.fasta

