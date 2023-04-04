# Amy Campbell
# concatenating and sorting forward & reverse reads

readpath="/home/acampbe/DFU/data/WGS_2020/SimulateGenomes141/SimReads"

for fwd in $readpath/*_R1.fastq ; do 
 read1="_R1.fastq"
 read2="_R2.fastq"
 nodirection=".fastq"
 sortedstring="sorted.fastq"

 rev=${fwd/$read1/$read2}
 combined=${fwd/$read1/$nodirection}
 combinedsorted=${fwd/$read1/$sortedstring}

 cat $fwd $rev > $combined
 cat $combined | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $combinedsorted



done
