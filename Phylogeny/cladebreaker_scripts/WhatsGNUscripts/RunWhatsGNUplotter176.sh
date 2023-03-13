# Amy Campbell
# Running case control plotter from WhatsGNU on:
# all DORNs except 1696 as cases
# DORN1696 + GCA_003336555.1, GCA_900097905.1, GCA_000248655.2, GCA_000149015.1 as controls
# to identify most differentiating CDS between these groups 

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

GNUReports=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/WG_reports/
CaseControl=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_176/case_controls_patient176.csv

WhatsGNU_plotter.py -st $CaseControl Patient176 $GNUReports
