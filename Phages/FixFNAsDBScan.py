from Bio import SeqIO
import os
import pandas
import sys
# Add any sequences that didn't get printed by DBSCAN-SWA despite
# being predicted as phage sequences (indexing bug in DBSCAN-SWA I think)
# also, for cases where it looks like a phage might have been interrupted by
# Unicycler's inference of an ORI, make a 'wraparound' FNA entry for that genome
# So that phage doesn't go undetected


# # Test file paths
# WrapPath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/TestingParse/WrapArounds.csv"
# MissingPath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/TestingParse/MissingSeqs.csv"
# OutputPath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/TestingParse/FNAsFixed"
# InputPath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/TestingParse/FNAs"
# ContigPath="/Users/amycampbell/Documents/DataInputGithub/data/IntraPatient/Phages/TestingParse/FinalSetIsolates"

WrapPath = sys.argv[1]
MissingPath = sys.argv[2]
OutputPath = sys.argv[3]
InputPath = sys.argv[4]
ContigPath = sys.argv[5]


filelist = os.listdir(InputPath)

WrapDF = pandas.read_csv(WrapPath)
WrapDF = WrapDF.sort_values(by="Start")

MissingDF = pandas.read_csv(MissingPath)
MissingDF = MissingDF[MissingDF['Missing']>0]
# will be true if DORN933 is among genomes that need wraps; False otherwise
########
#Contains933 = (WrapDF['Genome'].eq('DORN933').any())

# the following just gets the rows where genome is DORN933
#######################
# DORN933Wrap = (WrapDF[WrapDF['Genome']=='DORN933'])

# Then this gets the first row of those, and getse the value of 'Start' in that row
#############
#StartSecondFrag = (DORN933Wrap.iloc[0].Start)

copyscriptlines =[]

for filename in filelist:
#for filename in ["DORN1523_DBSCAN-SWA_prophage.fna"]:

    filepath = os.path.join(InputPath, filename)
    GenomeID = filename.replace("_DBSCAN-SWA_prophage.fna", "")
    newfilename=GenomeID + "_fixed.fna"
    Wrap = (WrapDF['Genome'].eq(GenomeID).any())
    MissingBool = (MissingDF['Genome'].eq(GenomeID).any())

    # ignore this genome if it doesn't need modification-- just set up a cp command for a shellscript for it
    if Wrap==False and MissingBool==False:
        linetoadd="cp " + str(os.path.join(InputPath, filename ) + " " + os.path.join(OutputPath, newfilename)) + "\n"
        copyscriptlines.append(linetoadd)

    else:
        SequenceDict = dict()

        records = list(SeqIO.parse(filepath, "fasta"))
        for i in range(len(records)):
            idstring = records[i].id
            sequence=records[i].seq
            idstringlist = (idstring.split('|'))
            contig=idstringlist[0]
            rangeseqlist = idstringlist[1].split(":")
            startseq = int(rangeseqlist[0])

            endseq = int(rangeseqlist[1])

            if(sequence==""):
                #print(startseq)
                #print(endseq)
                #print(contig)
                FullFasta = SeqIO.parse(os.path.join(ContigPath, (GenomeID + "_Final.fasta")), "fasta")
                contigfound=False
                while(contigfound==False):
                    contigrecord = next(FullFasta)
                    if contigrecord.id==contig:
                        NewSequence = (contigrecord.seq)[startseq:(endseq+1)]
                        SequenceDict[idstring]=NewSequence
                        contigfound=True

            else:
                SequenceDict[idstring]=sequence

        if Wrap==True:
            SubsetRows = WrapDF[WrapDF['Genome']==GenomeID]
            SecondFragmentEnd=int(SubsetRows.iloc[0].End)
            FirstFragmentStart=int(SubsetRows.iloc[1].Start)
            ContigToLook=(SubsetRows.iloc[0].Contig)
            FullFastaAgain = SeqIO.parse(os.path.join(ContigPath, (GenomeID + "_Final.fasta")), "fasta")
            sliceContig_found=False

            idstring_concat = str(ContigToLook) + "|" + str(FirstFragmentStart) + ":" + str(SecondFragmentEnd) + "|Concatenated"
            print(GenomeID)
            while(sliceContig_found==False):
                slicecontig_record = next(FullFastaAgain)
                if str(slicecontig_record.id) == str(ContigToLook):
                    firstfragment=str(slicecontig_record.seq[FirstFragmentStart:])
                    secondfragment=str(slicecontig_record.seq[:(SecondFragmentEnd+1)])
                    concat_fragment=firstfragment + secondfragment
                    print(idstring_concat)
                    SequenceDict[idstring_concat] = concat_fragment
                    sliceContig_found=True
        with open(os.path.join(OutputPath, newfilename), "w") as outputfile:
            for k in SequenceDict.keys():
                outputfile.write(k)
                outputfile.write('\n')
                outputfile.write(str(SequenceDict[k]))
                outputfile.write('\n')


ShellScriptCopy=open("CopyOKfnas.sh", "w")
ShellScriptCopy.writelines(copyscriptlines)
ShellScriptCopy.close()
