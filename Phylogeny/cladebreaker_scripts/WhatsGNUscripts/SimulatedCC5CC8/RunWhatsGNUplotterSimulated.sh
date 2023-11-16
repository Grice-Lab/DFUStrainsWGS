# Amy Campbell
# Running case control plotter from WhatsGNU on:
# CC8(case) and CC5(control) representative genomes
# to identify most differentiating CDS between these groups 

source ~/mambaforge/bin/activate ~/mambaforgeOLD/envs/WhatsGNUEnv

GNUReports=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/WG_reports/
CaseControl=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/case_controls_simulated.csv

WhatsGNU_plotter.py -st $CaseControl PatientSimulated $GNUReports
