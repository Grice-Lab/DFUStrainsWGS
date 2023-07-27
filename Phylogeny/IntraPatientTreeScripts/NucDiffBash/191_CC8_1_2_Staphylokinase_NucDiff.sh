#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2093_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2221_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2093_DORN2221 DORN2093_DORN2221

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2093_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1946_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2093_DORN1946 DORN2093_DORN1946

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2093_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2130_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2093_DORN2130 DORN2093_DORN2130

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2093_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2058_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2093_DORN2058 DORN2093_DORN2058

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2093_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1961_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2093_DORN1961 DORN2093_DORN1961

