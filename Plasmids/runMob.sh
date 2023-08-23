#!/bin/bash
# Amy Campbell
# Using MOB-SUITE to type the plasmids in the DFU100 
# 220 sequenced S. aureus isolates


source ~/mambaforge/bin/activate PlasmidEnv

#outputdir1214="/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/DORN1214"

plasmidsoutput=/home/acampbe/DFU/data/WGS_2020/MOB_Plasmids/
mkdir -p $plasmidsoutput

fastaext="_Final.fasta"
blank=""

for isolate in /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN*; do
	basenameisolate=$(basename $isolate)
	prefixisolate=${basenameisolate/$fastaext$blank}
	mob_recon -u --infile $isolate --outdir $plasmidsoutput$prefixisolate

done
