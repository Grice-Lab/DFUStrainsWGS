#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1399_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1410_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1399_DORN1410 DORN1399_DORN1410

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1399_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1413_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1399_DORN1413 DORN1399_DORN1413

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1399_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN1416_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN1399_DORN1416 DORN1399_DORN1416

