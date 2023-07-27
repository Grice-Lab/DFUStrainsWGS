#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN933_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN881_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN933_DORN881 DORN933_DORN881

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN933_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN880_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN933_DORN880 DORN933_DORN880

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN933_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1194_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN933_DORN1194 DORN933_DORN1194

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN933_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN976_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN933_DORN976 DORN933_DORN976

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN933_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN882_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN933_DORN882 DORN933_DORN882

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN933_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN925_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN933_DORN925 DORN933_DORN925

