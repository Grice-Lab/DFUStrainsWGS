# Amy Campbell
# 3/11/23 Running WhatsGNU main script on the reference genomes I want to include for patient 176 plotter call

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

StaphDB=/home/acampbe/DownloadedDatabases/WhatsGNU_Sau_Ortholog/Sau_Ortholog_10350.pickle
#roaryclusters=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/roary_alignment/results/clustered_proteins

#rm -r /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_003336555.1/WhatsGNU/
#output176_GCA_003336555=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_003336555.1/WhatsGNU/
#input176_GCA_003336555=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_003336555.1/annotation/GCA_003336555.1.faa

#rm -r /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_900097905.1/WhatsGNU/
#output176_GCA_900097905=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_900097905.1/WhatsGNU/
#input176_GCA_900097905=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_900097905.1/annotation/GCA_900097905.1.faa

#rm -r /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000248655.2/WhatsGNU/
#output176_GCA_000248655=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000248655.2/WhatsGNU/
#input176_GCA_000248655=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000248655.2/annotation/GCA_000248655.2.faa

#rm -r /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000149015.1/WhatsGNU/
#output176_GCA_000149015=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000149015.1/WhatsGNU/
#input176_GCA_000149015=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000149015.1/annotation/GCA_000149015.1.faa


#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_003336555 $input176_GCA_003336555
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_900097905 $input176_GCA_900097905
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_000248655 $input176_GCA_000248655
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_000149015 $input176_GCA_000149015


#GCA_900097905.1_WhatsGNU_report.txt
mkdir -p /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/WG_DORNs/
outputfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/WG_DORNs/
inputfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/

outputstring=DORN1602
inputstring=DORN1602/annotation/DORN1602.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1643
inputstring=DORN1643/annotation/DORN1643.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1644
inputstring=DORN1644/annotation/DORN1644.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1645
inputstring=DORN1645/annotation/DORN1645.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1646
inputstring=DORN1646/annotation/DORN1646.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1663
inputstring=DORN1663/annotation/DORN1663.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1695
inputstring=DORN1695/annotation/DORN1695.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1696
inputstring=DORN1696/annotation/DORN1696.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring

outputstring=DORN1731
inputstring=DORN1731/annotation/DORN1731.faa
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $outputfolder$outputstring $inputfolder$inputstring


