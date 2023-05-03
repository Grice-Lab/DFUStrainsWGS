







#!bin/bash
source /home/acampbe/mambaforge/bin/activate MetagenomicCladeEnv
Output=/home/acampbe/DFU/data/StressOperonMarkers/VariantComparisons

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_sigB.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_4_LAC_sigB.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_4.txt
Transcriptstring=SAUSA300_2022

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbW.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_4_LAC_rsbW.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_4.txt
Transcriptstring=SAUSA300_2023

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbV.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_4_LAC_rsbV.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_4.txt
Transcriptstring=SAUSA300_2024

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbU.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_4_LAC_rsbU.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_4.txt
Transcriptstring=SAUSA300_2025

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring


Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_sigB.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_6_LAC_sigB.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_6.txt
Transcriptstring=SAUSA300_2022

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbW.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_6_LAC_rsbW.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_6.txt
Transcriptstring=SAUSA300_2023

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbV.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_6_LAC_rsbV.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_6.txt
Transcriptstring=SAUSA300_2024

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbU.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_6_LAC_rsbU.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_6.txt
Transcriptstring=SAUSA300_2025

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring


Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_sigB.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_9_LAC_sigB.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_9.txt
Transcriptstring=SAUSA300_2022

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbW.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_9_LAC_rsbW.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_9.txt
Transcriptstring=SAUSA300_2023

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbV.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_9_LAC_rsbV.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_9.txt
Transcriptstring=SAUSA300_2024

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbU.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_9_LAC_rsbU.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_9.txt
Transcriptstring=SAUSA300_2025

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring


Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_sigB.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_12_LAC_sigB.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_12.txt
Transcriptstring=SAUSA300_2022

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbW.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_12_LAC_rsbW.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_12.txt
Transcriptstring=SAUSA300_2023

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbV.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_12_LAC_rsbV.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_12.txt
Transcriptstring=SAUSA300_2024

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

Early=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_1_LAC_rsbU.vcf
Late=/home/acampbe/DFU/data/StressOperonMarkers/AnnotatedSNPs/kalan01_DFUwgs_124_12_LAC_rsbU.vcf
BasesEarly=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_1.txt
BasesLate=/home/acampbe/DFU/data/StressOperonMarkers/Alignments/BaseCounts/kalan01_DFUwgs_124_12.txt
Transcriptstring=SAUSA300_2025

Rscript ../SnpEffOutput.R $Early $Late $BasesEarly $BasesLate $Output $Transcriptstring

