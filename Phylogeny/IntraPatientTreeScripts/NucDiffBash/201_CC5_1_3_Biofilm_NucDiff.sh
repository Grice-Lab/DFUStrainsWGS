#!bin/bash
# Making nucdiff scripts for intrapatient comparisons
mamba ~/mambaforge/bin/activate NucDiffEnv

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2213_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2105_DORN2213 DORN2105_DORN2213

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2137_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2105_DORN2137 DORN2105_DORN2137

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2148_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2105_DORN2148 DORN2105_DORN2148

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2149_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2105_DORN2149 DORN2105_DORN2149

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2105_Final.fasta /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/DORN2178_Final.fasta /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/NucDiffOutput/DORN2105_DORN2178 DORN2105_DORN2178

