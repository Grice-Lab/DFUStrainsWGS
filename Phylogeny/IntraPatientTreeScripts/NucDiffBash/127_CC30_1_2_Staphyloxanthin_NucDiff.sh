#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN613_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN620_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN613_DORN620 DORN613_DORN620

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN613_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN596_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN613_DORN596 DORN613_DORN596

