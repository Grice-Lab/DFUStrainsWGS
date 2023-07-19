# Amy Campbell 
# 2023
# Step 1 RNASeq on Two SA strains: FastQC

source /home/acampbe/software/miniconda3/bin/activate TrimmingEnvironment

READSPATH=/project/grice/storage/RNASeq_JE2_NE404/SeqCenterReads
OUTPUT=/project/grice/storage/RNASeq_JE2_NE404/FastQC
mkdir -p $OUTPUT

r1string="_R1.fastq.gz"
r2string="_R2.fastq.gz"
for filename in $READSPATH/*_R1.fastq.gz; do
        run1=$filename  
        run2=${run1/$r1string/$r2string}
        
        fastqc -o $OUTPUT $run1 -f fastq
        fastqc -o $OUTPUT $run2 -f fastq

done    

