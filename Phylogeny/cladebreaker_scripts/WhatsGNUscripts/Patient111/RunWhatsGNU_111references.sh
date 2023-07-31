# Amy Campbell
# 3/26/23 Running WhatsGNU main script on the reference genomes I want to include for patient 111 plotter call

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

StaphDB=/home/acampbe/DownloadedDatabases/WhatsGNU_Sau_Ortholog/Sau_Ortholog_10350.pickle
roaryclusters=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/roary_alignment/results/clustered_proteins
inputfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/
outputfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/WG_reports/

mkdir -p $outputfolder

# Move all WG reports generated by cladebreaker into the WG report folder, and then remove 1334 because we're ignoring it for this exercise
find $inputfolder -name '*WhatsGNU_report.txt' | xargs cp -t $outputfolder

# 1. GCA_900097905.1
input1=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_900097905.1/annotation/GCA_900097905.1.faa
output1=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_900097905.1/WhatsGNU
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output1 $input1
wgreport1="/GCA_900097905.1_WhatsGNU_report.txt"
cp $output1$wgreport1 $outputfolder


# 2. GCA_000248655.2 
input2=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_000248655.2/annotation/GCA_000248655.2.faa
output2=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_000248655.2/WhatsGNU
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output2 $input2
wgreport2="/GCA_000248655.2_WhatsGNU_report.txt"                
cp $output2$wgreport2 $outputfolder



# 3. GCA_000577355.1
input3=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_000577355.1/annotation/GCA_000577355.1.faa
output3=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_000577355.1/WhatsGNU
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output3 $input3
wgreport3="/GCA_000577355.1_WhatsGNU_report.txt"
cp $output3$wgreport3 $outputfolder


# 4. GCA_003336555.1
input4=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_003336555.1/annotation/GCA_003336555.1.faa
output4=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_003336555.1/WhatsGNU
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output4 $input4
wgreport4="/GCA_003336555.1_WhatsGNU_report.txt"
cp $output4$wgreport4 $outputfolder


# 5. GCA_900082315.1
input5=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_900082315.1/annotation/GCA_900082315.1.faa
output5=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_900082315.1/WhatsGNU
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output5 $input5
wgreport5="/GCA_900082315.1_WhatsGNU_report.txt"
cp $output5$wgreport5 $outputfolder

# 6. GCA_900081385.1
input6=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_900081385.1/annotation/GCA_900081385.1.faa
output6=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_900081385.1/WhatsGNU
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output6 $input6
wgreport6="/GCA_900081385.1_WhatsGNU_report.txt"
cp $output6$wgreport6 $outputfolder

# 7. GCA_003813255.1
input7=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_003813255.1/annotation/GCA_003813255.1.faa
output7=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/GCA_003813255.1/WhatsGNU
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output7 $input7
wgreport7="/GCA_003813255.1_WhatsGNU_report.txt"
cp $output7$wgreport7 $outputfolder


