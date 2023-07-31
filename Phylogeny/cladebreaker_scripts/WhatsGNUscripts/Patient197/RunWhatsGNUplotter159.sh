# Amy Campbell
# Running case control plotter from WhatsGNU on:
# CC15 ‘Case’: DORN1364, DORN1384, DORN1340, DORN1339, DORN1352
# CC1 ‘Control’: DORN1353,GCA_000149015.1,GCA_000248655.2 , GCA_003336555.1, GCA_000577355.1
# to identify most differentiating CDS between these groups in patient 159

source ~/mambaforge/bin/activate ~/mambaforge/envs/WhatsGNUEnv

GNUReports=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/WG_reports/
CaseControl=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/case_controls_patient159.csv

WhatsGNU_plotter.py -st $CaseControl Patient159 $GNUReports
