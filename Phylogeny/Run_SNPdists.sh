#!/bin/bash
# get SNP distances from the core genome alignments output by Roary

source /home/acampbe/software/miniconda3/bin/activate prokenv

snp-dists /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/core_gene_alignment.aln > /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/snp-dists.tsv

snp-dists /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/core_gene_alignment.aln > /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput_withRefs/snp-dists-References.tsv
