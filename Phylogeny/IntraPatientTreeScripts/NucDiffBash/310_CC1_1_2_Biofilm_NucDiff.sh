#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1082_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1085_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1082_DORN1085 DORN1082_DORN1085

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1082_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1081_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1082_DORN1081 DORN1082_DORN1081

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1082_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1086_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1082_DORN1086 DORN1082_DORN1086

