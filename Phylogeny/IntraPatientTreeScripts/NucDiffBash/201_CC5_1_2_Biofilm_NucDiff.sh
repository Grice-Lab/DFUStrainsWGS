#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2148_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2148 DORN2136_DORN2148

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2105_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2105 DORN2136_DORN2105

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2149_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2149 DORN2136_DORN2149

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2178_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2178 DORN2136_DORN2178

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2179_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2179 DORN2136_DORN2179

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2136_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2150_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2136_DORN2150 DORN2136_DORN2150

