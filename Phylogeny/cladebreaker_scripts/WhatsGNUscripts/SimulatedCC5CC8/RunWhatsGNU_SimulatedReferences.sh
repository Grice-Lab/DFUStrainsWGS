# Amy Campbell
# 3/11/23 Running WhatsGNU main script on the reference genomes I want to include for patient 176 plotter call

source ~/mambaforge/bin/activate ~/mambaforgeOLD/envs/WhatsGNUEnv

StaphDB=/home/acampbe/DownloadedDatabases/WhatsGNU_Sau_Ortholog/Sau_Ortholog_10350.pickle


outerfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/
wgreportsfolder=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/WG_reports/

mkdir -p $wgreportsfolder

cp /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/SaureusCC5/WhatsGNU/WhatsGNU_Report/SaureusCC5_WhatsGNU_report.txt $wgreportsfolder
cp /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/SaureusCC8/WhatsGNU/WhatsGNU_Report/SaureusCC8_WhatsGNU_report.txt	$wgreportsfolder

output1=$outerfolder"GCA_000010465.1/WhatsGNU/"
input1=$outerfolder"GCA_000010465.1/annotation/GCA_000010465.1.faa"

output2=$outerfolder"GCA_002310395.1/WhatsGNU/"
input2=$outerfolder"GCA_002310395.1/annotation/GCA_002310395.1.faa"

output3=$outerfolder"GCA_900092595.1/WhatsGNU/"
input3=$outerfolder"GCA_900092595.1/annotation/GCA_900092595.1.faa"

output4=$outerfolder"GCA_900457575.1/WhatsGNU/"
input4=$outerfolder"GCA_900457575.1/annotation/GCA_900457575.1.faa"

output5=$outerfolder"GCA_001019305.2/WhatsGNU/"
input5=$outerfolder"GCA_001019305.2/annotation/GCA_001019305.2.faa"


output6=$outerfolder"GCA_001611465.1/WhatsGNU/"
input6=$outerfolder"GCA_001611465.1/annotation/GCA_001611465.1.faa"


output7=$outerfolder"GCA_001611475.1/WhatsGNU/"
input7=$outerfolder"GCA_001611475.1/annotation/GCA_001611475.1.faa"


output8=$outerfolder"GCA_001611545.1/WhatsGNU/"
input8=$outerfolder"GCA_001611545.1/annotation/GCA_001611545.1.faa"


#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output1 $input1
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output2 $input2
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output3 $input3
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output4 $input4


#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output5 $input5
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output6 $input6
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output7 $input7
#WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output8 $input8

find $outerfolder -name \*_WhatsGNU_report.txt -exec cp {} /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/WG_reports/ \;
