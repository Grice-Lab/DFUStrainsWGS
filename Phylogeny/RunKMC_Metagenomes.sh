#!bin/bash
# Amy Campbell
# Testing for presence/absence of CCs in metagenomic samples 
# By checking for CC-unique 41-mers

source ~/mambaforge/bin/activate ~/mambaforge/envs/pankmer

# variables that will be positional in final form 
InputFastqPath=/home/acampbe/DFU/data/DFU_Metagenome_Microbes/
OutputFolder=/home/acampbe/DFU/data/CladeRepresentatives/MetagenomeKmers/
CCkmersPath="/home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/KmerOutput/"
ResultsOutput="/home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/SampleKmerCounts/"
TsvOutput="/home/acampbe/DFU/data/WGS_2020/CladeRepresentatives/TSVs/"


mkdir -p $OutputFolder
mkdir -p $ResultsOutput
mkdir -p $TsvOutput


# some string variables
prefixstring="_Unique.kmc_pre"
blankstring=""
Fastqext=".fastq.gz"
kmersext="_kmers"
underscore="_"

#fastqfile="/home/acampbe/DFU/data/DFU_Metagenome_Microbes/kalan01_DFUwgs_124_6.fastq.gz"

for fastqfile in $InputFastqPath* ; do 

	fastqfileBase=$(basename $fastqfile)
	fastqID=${fastqfileBase/$Fastqext/$blankstring}


	KmerOutputMetagenome=$OutputFolder$fastqID$kmersext

	kmc -k41 -fq $fastqfile $KmerOutputMetagenome .

	for CCfile in $CCkmersPath*_Unique.kmc_pre; do
    		# Get the pre and suf db files for the CC of interest 
    		prefile=$CCfile
    		suffile=${CCfile/pre/suf}
    		preext=".kmc_pre"

    		NoExt=${CCfile/$preext/$blankstring}
    		echo $NoExt
    		basestring=$(basename $CCfile)
    		CCid=${basestring/$prefixstring/$blankstring}
    		outputfile=$ResultsOutput$fastqID$underscore$CCid
    		outputtxt=$outputfile".txt"
    		# Call kmc_tools to find intersection of 
    		kmc_tools simple $NoExt $KmerOutputMetagenome intersect $outputfile -ocright
    
    		kmc_tools transform $outputfile dump $outputtxt
    
	done

done
rm $KmerOutputMetagenome*


exttsv="_KMerComposition.tsv"
SampleTSV=$TsvOutput$fastqID$exttsv
textext=".txt"
tab="\t"
tempext="_temp.txt"
tempext2="_temp2.txt"

touch $SampleTSV
for textfile in $ResultsOutput$fastqID$underscore*.txt; do 

    baseText=$(basename $textfile)
    CC_string=${baseText/$fastqID$underscore/$blankstring}
    CC_string=${CC_string/$textext/$blankstring}
    StringPrefix=$CC_string$tab
    tempfile=${textfile/$textext/$tempext}
    tempfile2=${textfile/$textext/$tempext2}
    awk -F '\t' '{ print $2 }' $textfile > $tempfile 
    awk -v prefix=$StringPrefix '{print prefix $0}' $tempfile > $tempfile2

    cat $tempfile2 >> $SampleTSV

    #rm $tempfile2
    #rm $tempfile

done



#/home/acampbe/DFU/data/CladeRepresentatives/SampleKmerCounts/kalan01_DFUwgs_124_6_CC7.txt
