# Amy Campbell
# Running case control plotter from WhatsGNU on CC30 and CC1 genomes
# to identify most differentiating CDS between these groups in patient 145

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

GNUReports=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/WG_reports/
CaseControl=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_145/case_controls_patient145.csv

WhatsGNU_plotter.py -st $CaseControl Patient145 $GNUReports
