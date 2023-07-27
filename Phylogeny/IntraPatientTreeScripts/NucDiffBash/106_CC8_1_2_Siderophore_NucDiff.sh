#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN76_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN105_DORN76 DORN105_DORN76

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN62_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN105_DORN62 DORN105_DORN62

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN47_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN105_DORN47 DORN105_DORN47

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN56_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN105_DORN56 DORN105_DORN56

