#!/bin/bash
# 
mkdir -p /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/
source /home/acampbe/software/miniconda3/bin/activate prokenv

# running Prokka annotation on the reference genomes for the Staph aureus tree with DORN samples 

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC1_MSSA476 --force --prefix CC1_MSSA476 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC1_MSSA476.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC1_MW2 --force --prefix CC1_MW2 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC1_MW2.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC22_HO --force --prefix CC22_HO --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC22_HO.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC30_MRSA252 --force --prefix CC30_MRSA252 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC30_MRSA252.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC398_ATCC6538 --force --prefix CC398_ATCC6538 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC398_ATCC6538.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC398 --force --prefix CC398 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC398.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC5_Mu50 --force --prefix CC5_Mu50 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC5_Mu50.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC5_N315 --force --prefix CC5_N315 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC5_N315.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC72_CN1 --force --prefix CC72_CN1 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC72_CN1.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC8_NCTC8325 --force --prefix CC8_NCTC8325 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC8_NCTC8325.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/CC8_Newman --force --prefix CC8_Newman --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/CC8_Newman.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/S_epidermidis --force --prefix S_epidermidis --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/S_epidermidis.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/USA100_AR465 --force --prefix USA100_AR465 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/USA100_AR465.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/USA300_FPR3757 --force --prefix USA300_FPR3757 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/USA300_FPR3757.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/USA400_051 --force --prefix USA400_051 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/USA400_051.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_502A --force --prefix SA_502A --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/SA_502A.fasta


prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_AR464 --force --prefix SA_AR464 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/newrefs/SA_AR464.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_CFSAN007883 --force --prefix SA_CFSAN007883 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/newrefs/SA_CFSAN007883.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_MCRF184 --force --prefix SA_MCRF184 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/newrefs/SA_MCRF184.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_AR464 --force --prefix SA_AR464 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/newrefs/SA_AR464.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_CFSAN007883 --force --prefix SA_CFSAN007883 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/newrefs/SA_CFSAN007883.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_MCRF184 --force --prefix SA_MCRF184 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/newrefs/SA_MCRF184.fasta

prokka --outdir /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/SA_UP_1150 --force --prefix SA_UP_1150 --genus Staphylococcus /home/acampbe/DFU/data/WGS_2020/RoaryResults/SAReferences/SA_UP_1150.fasta


