# Amy Campbell
# Running case control plotter from WhatsGNU on:
# DORN2083, DORN2144, DORN2089, DORN2123, DORN2187  as control
# DORN2205, GCA_000606925.1, GCA_900082865.1, GCA_001063675.1, GCA_001065215.1, DORN2205 as case
# to identify most differentiating CDS between these groups 

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

GNUReports=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/WG_reports/

find /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/ -name '*WhatsGNU_report.txt' | xargs cp -t /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/WG_reports/

#/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/WG_reports/

#rm -r /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/Patient197
CaseControl=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_197/case_controls_patient197.csv

WhatsGNU_plotter.py -st $CaseControl Patient197 $GNUReports
