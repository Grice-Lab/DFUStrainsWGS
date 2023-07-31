# Amy Campbell
# Running case control plotter from WhatsGNU on CC30 and CC1 genomes
# to identify most differentiating CDS between these groups in patient 111

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

GNUReports=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/WG_reports/
CaseControl=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_111/case_controls_patient111.csv

WhatsGNU_plotter.py -st $CaseControl Patient111 $GNUReports
