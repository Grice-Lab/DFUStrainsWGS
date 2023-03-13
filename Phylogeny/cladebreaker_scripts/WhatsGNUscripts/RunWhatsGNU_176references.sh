# Amy Campbell
# 3/11/23 Running WhatsGNU main script on the reference genomes I want to include for patient 176 plotter call

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

StaphDB=/home/acampbe/DownloadedDatabases/WhatsGNU_Sau_Ortholog/Sau_Ortholog_10350.pickle

output176_GCA_003336555=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_003336555.1/WhatsGNU/
input176_GCA_003336555=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_003336555.1/annotation/GCA_003336555.1.faa

output176_GCA_900097905=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_900097905.1/WhatsGNU/
input176_GCA_900097905=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_900097905.1/annotation/GCA_900097905.1.faa

output176_GCA_000248655=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000248655.2/WhatsGNU/
input176_GCA_000248655=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000248655.2/annotation/GCA_000248655.2.faa

output176_GCA_000149015=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000149015.1/WhatsGNU/
input176_GCA_000149015=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/GCA_000149015.1/annotation/GCA_000149015.1.faa



WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_003336555 $input176_GCA_003336555
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_900097905 $input176_GCA_900097905
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_000248655 $input176_GCA_000248655
WhatsGNU_main.py -d $StaphDB -dm ortholog -o $output176_GCA_000149015 $input176_GCA_000149015



