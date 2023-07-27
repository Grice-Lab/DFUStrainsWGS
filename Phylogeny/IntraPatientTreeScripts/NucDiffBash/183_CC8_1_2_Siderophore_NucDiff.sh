#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1743_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1782_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1743_DORN1782 DORN1743_DORN1782

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1743_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1729_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1743_DORN1729 DORN1743_DORN1729

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1743_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1747_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1743_DORN1747 DORN1743_DORN1747

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1743_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1765_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1743_DORN1765 DORN1743_DORN1765

