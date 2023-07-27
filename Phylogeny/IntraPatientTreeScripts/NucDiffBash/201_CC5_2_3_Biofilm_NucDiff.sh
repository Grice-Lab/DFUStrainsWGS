#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2213_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2213 DORN2136_DORN2213

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2137_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2137 DORN2136_DORN2137

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2179_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2179 DORN2136_DORN2179

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2150_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2150 DORN2136_DORN2150

