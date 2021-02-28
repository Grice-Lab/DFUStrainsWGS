#!/bin/bash
# Amy Campbell
# Feb 2021
# Run eggnog emapper to annotate roary-output pan_genome_reference.fa sequences 
# with eggnog database (including COG functional categories

source /home/acampbe/software/miniconda3/bin/activate eggNOGenv
export EGGNOG_DATA_DIR=/home/acampbe/EggNOGdb/data
export PATH=/home/acampbe/software/eggnog-mapper-2.0.8/eggnogmapper/bin:"$PATH"
python3 emapper.py -m diamond -i /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/pan_genome_reference.fa --itype CDS -o panGenome207
mv *panGenome207* /home/acampbe/DFU/data/WGS_2020/EggnogOutput/
