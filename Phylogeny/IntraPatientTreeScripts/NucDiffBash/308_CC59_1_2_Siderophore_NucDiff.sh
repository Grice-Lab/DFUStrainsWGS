#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN968_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN962_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN968_DORN962 DORN968_DORN962

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN968_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN915_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN968_DORN915 DORN968_DORN915

