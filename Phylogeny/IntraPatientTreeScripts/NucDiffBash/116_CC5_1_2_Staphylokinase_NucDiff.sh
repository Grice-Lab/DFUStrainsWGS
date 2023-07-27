#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN283_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN280_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN283_DORN280 DORN283_DORN280

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN283_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN300_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN283_DORN300 DORN283_DORN300

