#!bin/bash
# Amy Campbell
# March 2023
# Take the output of the alignments of metagenomic reads to patient-specific markers
# calls CalculateComposition.R to calculate the composition of that patient's wound S. aureus over time

source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv


# Patient 176
##############
BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient176/output/BCFs"
VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/Patient176/Patient176_Variants.tsv
Patient="176"
Case="CC5"
Control="CC1"
OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions"

#Rscript CalculateComposition.R $BCFs $VariantListPath $Patient $Case $Control $OutputFolder

# Patient 159
##############
#BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient159/output/BCFs"
#VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/Patient159_Variants.tsv
Patient="159"
Case="CC15"
Control="CC1"
OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions"
#Rscript CalculateComposition.R $BCFs $VariantListPath $Patient $Case $Control $OutputFolder

# Patient 159
##############


BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient145/output/BCFs/"
mkdir -p $BCFs
#mv /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient145/output/*.bcf $BCFs
VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/Patient145/Patient145_Variants.tsv
Patient="145"
Case="CC30"
Control="CC1"
OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions"
Rscript CalculateComposition.R $BCFs $VariantListPath $Patient $Case $Control $OutputFolder


# Patient 197
##############

BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient197/output/BCFs/"
mkdir -p $BCFs
#mv /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient197/output/*.bcf $BCFs
VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197/Patient197_Variants.tsv
Patient="197"
Case="CC8_1"
Control="CC8_2"
OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions"
#Rscript CalculateComposition.R $BCFs $VariantListPath $Patient $Case $Control $OutputFolder

# Patient 191
##############

BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/output/BCFs/"
mkdir -p $BCFs
#mv /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient191/output/*.bcf $BCFs
VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_191/Patient191/Patient191_Variants.tsv
Patient="191"
Case="CC5"
Control="CC8"
OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions"
#Rscript CalculateComposition.R $BCFs $VariantListPath $Patient $Case $Control $OutputFolder

# Patient 173
##############

#BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient173/output/BCFs/"
#mkdir -p $BCFs
#mv /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient173/output/*.bcf $BCFs
#VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_173/Patient173/Patient173_Variants.tsv
#Patient="173"
#Case="CC72"
#Control="CC5"
#OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions"
#Rscript CalculateComposition.R $BCFs $VariantListPath $Patient $Case $Control $OutputFolder


# Patient 111
##############

BCFs="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/output/BCFs/"
#mkdir -p $BCFs
#mv /home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Patient111/output/*.bcf $BCFs
VariantListPath=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/Patient111/Patient111_Variants.tsv
Patient="111"
Case="CC1"
Control="CC5"
OutputFolder="/home/acampbe/DFU/data/AlignMetagenomesCladeMarkers/Compositions"
Rscript CalculateComposition.R $BCFs $VariantListPath $Patient $Case $Control $OutputFolder

