#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN468_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN460_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN468_DORN460 DORN468_DORN460

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN468_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN467_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN468_DORN467 DORN468_DORN467

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN468_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN459_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN468_DORN459 DORN468_DORN459

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN468_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN489_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN468_DORN489 DORN468_DORN489

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN468_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN498_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN468_DORN498 DORN468_DORN498

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN468_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN564_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN468_DORN564 DORN468_DORN564

