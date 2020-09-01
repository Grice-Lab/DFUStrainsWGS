#!/bin/bash
source /home/acampbe/software/miniconda3/bin/activate prokenv
roary -e -p -8 -f /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput /home/acampbe/DFU/data/WGS_2020/RoaryResults/gff_files/*
