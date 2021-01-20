# Amy Campbell
# December 2020
# Association testing with DBGWAS 

source /home/acampbe/software/miniconda3/bin/activate DBGWASEnv

# Original run(219 isolates)
# /home/acampbe/software/DBGWAS/bin/DBGWAS -strains XanthinMapping -newick /home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_bestTree.RaxMLTreeNoRefs.newick -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_WithAnnotations/


# k-mer size 31 (All 221 Isolates)
# /home/acampbe/software/DBGWAS/bin/DBGWAS -strains XanthinMappingFull -maf .05 -newick /home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_bestTree.RaxMLTreeNoRefs.newick -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_31 -verbose 3

# k-mer size 41
/home/acampbe/software/DBGWAS/bin/DBGWAS -strains XanthinMappingFull -k 41 -maf .05 -newick  /home/acampbe/DFU/data/WGS_2020/Trees/RaxML_NoRefs/RAxML_bestTree.RaxMLTreeNoRefs.newick -nc-db /home/acampbe/DFU/data/WGS_2020/RoaryResults/roaryoutput/pan_genome_reference.fa -output /home/acampbe/DFU/data/WGS_2020/DBGWAS_Output_Full_41

