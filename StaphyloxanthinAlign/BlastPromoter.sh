#!/bin/bash
# Blast searching for promoter for crtOPQMN operon in each of the assemblies
source /home/acampbe/software/miniconda3/bin/activate BlastEnv

#makeblastdb -in XanthinPromoter_BasedOnPrimers.fasta -input_type fasta -dbtype nucl -out XanthinPromoter

for fastaname in /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_Cleaned_Assemblies/FinalContigs/*.fasta; do
        base=$(basename $fastaname)
	baseext=".fasta"
	blank=""
	extensionless=${base/$baseext/$blank}
	tabext=".tab"
	makeblastdb -in $fastaname -dbtype nucl -out $extensionless
        blastn -query xanthinpromoter.fasta -out $extensionless$tabext -db $extensionless -outfmt "6 qseqid sseqid pident length sstart send qcovs evalue bitscore"
done
