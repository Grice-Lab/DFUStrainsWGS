#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1540_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1531_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1540_DORN1531 DORN1540_DORN1531

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1540_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1623_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1540_DORN1623 DORN1540_DORN1623

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1540_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1562_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1540_DORN1562 DORN1540_DORN1562

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1540_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1546_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1540_DORN1546 DORN1540_DORN1546

