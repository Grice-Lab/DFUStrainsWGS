#!/bin/bash 
# Amy Campbell
# 11-2020
# blastn searches of each circularized contig in the CC5 lineages against "plsdb" downloaded 2020-11-25
mkdir -p /home/acampbe/CC5Plasmids/BlastResults

#PLSDB citation: 
################
# Valentina Galata, Tobias Fehlmann, Christina Backes, Andreas Keller; 
#PLSDB: a resource of complete bacterial plasmids, Nucleic Acids Res., 2018 Oct 31, doi: 10.1093/nar/gky1050
export BLASTDB='/home/acampbe/DownloadedDatabases/BlastDBs'

source /home/acampbe/software/miniconda3/bin/activate BlastEnv
#for file in /home/acampbe/CC5Plasmids/*.fasta ; do
#	otputfolder="/home/acampbe/CC5Plasmids/BlastResults/"
#        otputstring=".out"
#        filebase=$(basename $file)
#        circext="_circular.fasta"
#        blank=""
#        noext=${filebase/$circext/$blank}

#        blastn -query $file -db plsdb.fna -out $otputfolder$noext$otptstring -outfmt "6 qseqid qlen slen pident evalue bitscore stitle"
#-outfmt "6 qseqid sseqid pident qcovs evalue bitscore"

#done

# CC8 
mkdir -p /home/acampbe/CC8Plasmids/BlastResults/
for file in /home/acampbe/CC8Plasmids/*.fasta ; do
        otputfolder="/home/acampbe/CC8Plasmids/BlastResults/"
        otputstring=".out"
        filebase=$(basename $file)
        circext="_circular.fasta"
        blank=""
        noext=${filebase/$circext/$blank}
	blastn -query $file -db plsdb.fna -out $otputfolder$noext$otptstring -outfmt "6 qseqid qlen slen pident evalue bitscore stitle"
        #blastn -query $file -db plsdb.fna -out $otputfolder$noext$otptstring -outfmt "6 qseqid qlen slen pident evalue bitscore $
#-outfmt "6 qseqid sseqid pident qcovs evalue bitscore"

done


