# Amy Campbell 
# 2023
# Trimming RNAseq reads for NE404 vs. JE2

source /home/acampbe/software/miniconda3/bin/activate TrimmingEnvironment

ReadsPath=/project/grice/storage/RNASeq_JE2_NE404/SeqCenterReads/
Output=/project/grice/storage/RNASeq_JE2_NE404/Trimmed
mkdir -p $Output

for filename in $ReadsPath/*_R1.fastq.gz; do
	r1string="_R1.fastq.gz"
	r2string="_R2.fastq.gz"
	blankstring=""

        run1=$filename 
        run2=${run1/$r1string/$r2string}
        basefname=$(basename $filename)
	basefnameNoExt=${basefname/$r1string/$blankstring}

        trim_galore --paired --clip_R1 1 --clip_R2 1 --basename $basefnameNoExt --output_dir $Output $run1 $run2

done    

