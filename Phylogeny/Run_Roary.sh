#!/bin/bash
# Updated May 2021

source /home/acampbe/software/miniconda3/bin/activate prokenv
roary -e -z -p -16 -f /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput /home/acampbe/DFU/data/WGS_2020/RoaryResults/gffs_wRefs_ContaminatedRemoved/*

