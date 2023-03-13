# Amy Campbell
# 3/11/23 Running WhatsGNU main script on the reference genomes I want to include for patient 197 plotter call

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

StaphDB=/home/acampbe/DownloadedDatabases/WhatsGNU_Sau_Ortholog/Sau_Ortholog_10350.pickle

output197_GCA_000606925=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_000606925.1/WhatsGNU/
input197_GCA_000606925=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_000606925.1/annotation/GCA_000606925.1.faa

output197_GCA_900082865=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_900082865.1/WhatsGNU/
input197_GCA_900082865=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_900082865.1/annotation/GCA_900082865.1.faa

output197_GCA_001063675=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_001063675.1/WhatsGNU/
input197_GCA_001063675=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_001063675.1/annotation/GCA_001063675.1.faa

output197_GCA_001065215=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_001065215.1/WhatsGNU/
input197_GCA_001065215=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/GCA_001065215.1/annotation/GCA_001065215.1.faa

WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output197_GCA_000606925 $input197_GCA_000606925
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output197_GCA_001063675 $input197_GCA_001063675
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output197_GCA_900082865 $input197_GCA_900082865
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output197_GCA_001065215 $input197_GCA_001065215

