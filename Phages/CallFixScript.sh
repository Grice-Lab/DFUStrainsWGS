#!/bin/bash
# Amy Campbell

source ~/mambaforge/bin/activate DBSCANSWA

outputfolder=/home/acampbe/DFU/data/WGS_2020/PhageResults/FixedFNAs
mkdir -p $outputfolder

inputfolder=/home/acampbe/DFU/data/WGS_2020/PhageResults/FNAs

contigpath=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates

wrappath=/home/acampbe/DFUStrainsWGS/Phages/WrapArounds.csv

missingpath=/home/acampbe/DFUStrainsWGS/Phages/MissingSeqs.csv

python3 FixFNAsDBScan.py $wrappath $missingpath $outputfolder $inputfolder $contigpath

