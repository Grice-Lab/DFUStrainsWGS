# Amy Campbell
# Running case control plotter from WhatsGNU for CC5 vs CC8, CC5 vs CC9
# to identify most differentiating CDS between these groups in patient 124

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

GNUReports=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/WG_reports/

CaseControl_cc9=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/case_controls_patient124_cc9.csv
CaseControl_cc8=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/case_controls_patient124_cc8.csv

WhatsGNU_plotter.py -st $CaseControl_cc8 Patient124_CC8 $GNUReports
WhatsGNU_plotter.py -st $CaseControl_cc9 Patient124_CC9 $GNUReports

# rm -r /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC8
# rm -r /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_124/Patient124_CC9
