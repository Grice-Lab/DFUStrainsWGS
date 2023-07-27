#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN56_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN76_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN56_DORN76 DORN56_DORN76

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN56_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN62_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN56_DORN62 DORN56_DORN62

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN56_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN47_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN56_DORN47 DORN56_DORN47

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN56_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN105_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN56_DORN105 DORN56_DORN105

