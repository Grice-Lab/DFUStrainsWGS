#!/bin/bash
# Amy campbell
# reference-based assembly of DORN1499, DORN1743, DORN1761, DORN1952, DORN2221
# Because it seems yabj/spoVG stretch was lost in their de novo assemblies


source ~/mambaforge/bin/activate YabJSpoVGEnv

outputfolder=/home/acampbe/DFU/data/WGS_2020/yabJ_spoVG/For_Reassembly/

# DORN1499
###########
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1499trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1499trimmedgalore_val_2.fastq
trustedcontigs=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1502_Final.fasta
outputfname=DORN1499reassembled.fasta
outputdir=DORN1499
fastaext=".fasta"
contigsfast="/contigs.fasta"
outputfname=$outputfname
spades.py -1 $readspath1 -2 $readspath2 --trusted-contigs $trustedcontigs -o $outputfolder$outputdir
cp $outputfolder$outputdir$contigsfast $outputfolder$outputdir/$outputfname

# DORN1743
###########
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1743trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1743trimmedgalore_val_2.fastq
trustedcontigs=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1782_Final.fasta
outputfname=DORN1743reassembled.fasta
outputdir=DORN1743
fastaext=".fasta"
contigsfast="/contigs.fasta"
outputfname=$outputfname        
spades.py -1 $readspath1 -2 $readspath2 --trusted-contigs $trustedcontigs -o $outputfolder$outputdir
cp $outputfolder$outputdir$contigsfast $outputfolder$outputdir/$outputfname 


# DORN1761
###########
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1761trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1761trimmedgalore_val_2.fastq
trustedcontigs=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1740_Final.fasta
outputfname=DORN1761reassembled.fasta
outputdir=DORN1761
fastaext=".fasta"
contigsfast="/contigs.fasta"
outputfname=$outputfname        
spades.py -1 $readspath1 -2 $readspath2 --trusted-contigs $trustedcontigs -o $outputfolder$outputdir
cp $outputfolder$outputdir$contigsfast $outputfolder$outputdir/$outputfname 


# DORN1952
###########
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1952trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN1952trimmedgalore_val_2.fastq
trustedcontigs=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2166_Final.fasta
outputfname=DORN1952reassembled.fasta
outputdir=DORN1952
fastaext=".fasta"
contigsfast="/contigs.fasta"
outputfname=$outputfname
spades.py -1 $readspath1 -2 $readspath2 --trusted-contigs $trustedcontigs -o $outputfolder$outputdir
cp $outputfolder$outputdir$contigsfast $outputfolder$outputdir/$outputfname

# DORN2221
###########
readspath1=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2221trimmedgalore_val_1.fastq
readspath2=/project/grice/storage/DFUShortReads2022/trimmedreads/DORN2221trimmedgalore_val_2.fastq
trustedcontigs=/home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2093_Final.fasta
outputfname=DORN2221reassembled.fasta
outputdir=DORN2221
fastaext=".fasta"
contigsfast="/contigs.fasta"
outputfname=$outputfname
spades.py -1 $readspath1 -2 $readspath2 --trusted-contigs $trustedcontigs -o $outputfolder$outputdir
cp $outputfolder$outputdir$contigsfast $outputfolder$outputdir/$outputfname

