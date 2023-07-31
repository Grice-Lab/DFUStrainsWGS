# Amy Campbell
# Takes in maps of core genes for each intra-patient comparison to their pangenomeIDs
# Makes a multifasta of each comparator genome's version of that core gene, which notes in id whether
# the genome is in the 'high' or 'low' cluster for the phenotype comparison

# also makes a call to do a mafft pairwise multiple alignment to the longest 'high' gene representative
# makes snp-sites call for the mafft output alignment

import os
from Bio import SeqIO
import pandas
import sys

# Takes in
PanGenomePath = sys.argv[1]  # path to roary intermediates' pan_genome_sequences
# e.g. /home/acampbe/DFU/data/WGS_2020/RoaryResultsPGAP2022/pan_genome_sequences/

CoreGeneMapPath=sys.argv[2] # path to all the core gene maps
# (e.g. /home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/CoreGeneMaps/)

OutputPathFastas =sys.argv[3]  # path to output fastas (/home/acampbe/DFU/data/WGS_2020/IntraPatientRoary/CoreGeneMultifastas)

OutputMAFFT_SNPsPath= sys.argv[4] # path to output snp-sites and mafft call files
# e.g., /home/acampbe/DFUStrainsWGS/Phylogeny/IntraPatientTreeScripts/MAFFT_SNPSites_Scripts

OutputMAFFT_SNPs_Outputs = sys.argv[5]
# e.g., /home/acampbe/DFUStrainsWGS/IntraPatientTreeScripts/MAFFT_SNPSites_Scripts


CoreGeneMapList = os.listdir(CoreGeneMapPath)

for mapitem in CoreGeneMapList:
    CoreGeneMap = os.path.join(CoreGeneMapPath,mapitem )
    comparisonsMAFFTcommands = []
    comparisonsSEDcommands = []
    comparisonsSNPcommands = []


    CoreGeneMapFileString = os.path.basename(CoreGeneMap).replace("CoreGeneMappings.csv", "")
    CoreGeneMapFileList = CoreGeneMapFileString.split("_")

    patientNum = (CoreGeneMapFileList[1])
    CC = (CoreGeneMapFileList[2])
    Cluster1 = (CoreGeneMapFileList[3])
    Cluster2 = (CoreGeneMapFileList[4])
    Phenotype=(CoreGeneMapFileList[5])

    CoreGeneDF = pandas.read_csv(CoreGeneMap)
    CoreGeneDF.columns = ['Ind','Gene','IsolateID','Phenotype', 'HighOrLow', 'PanGenomeID']



    numUniquegenes=0
    for gene in set(CoreGeneDF['Gene']):
        littleDF = CoreGeneDF[CoreGeneDF['Gene']==gene]
        geneReformatted = gene.replace('\'','_')
        geneReformatted = geneReformatted.replace('(','_')
        geneReformatted = geneReformatted.replace(')','_')
        geneReformatted = geneReformatted.replace(':','_')
        geneReformatted = geneReformatted.replace('-','_')

        genefastapath=os.path.join(PanGenomePath, str(geneReformatted) + ".fa.aln")
        
        sequencesgene = SeqIO.to_dict(SeqIO.parse(genefastapath, "fasta"))
        sequencesgene = {key:val for key, val in sequencesgene.items() if key in list(littleDF.PanGenomeID)}

        sequences = list(sequencesgene.values())
        sequences = list(map(lambda x: str(x.seq), sequences))

        outputfastapath = os.path.join(OutputPathFastas, geneReformatted + "_" + patientNum + "_" + CC + "_" + Cluster1 + "_" + Cluster2 + "_" + Phenotype  +".fasta")
    if (len(set(sequences))) > 1:

        numUniquegenes=numUniquegenes+1
        with open(outputfastapath, "w") as fastoutput:
            maxlength = 0
            for row in range(littleDF.shape[0]):
                rowitem = (littleDF.iloc[row,])
                keytouse = rowitem['PanGenomeID']
                genename = rowitem['Gene']
                high_low = rowitem['HighOrLow']
                isolate = rowitem['IsolateID']
                sequencestring = str(sequencesgene[keytouse].seq)
                sequencestring = sequencestring.replace("-", "")
                fastoutput.write(str(">"+isolate+"_" + high_low + "\n"))
                fastoutput.write(sequencestring + "\n")

                # make a reference fasta with the longest 'high'
                if (len(sequencestring) > maxlength) & (high_low=="High"):
                    referencegene = isolate
                    referenceseq =  sequencestring
                    maxlength = len(sequencestring)

        outputreferencefasta = os.path.join(OutputPathFastas,"Reference_"+ geneReformatted + "_" + patientNum + "_" + CC + "_" + Cluster1 + "_" + Cluster2 + "_" + Phenotype  +".fasta")
        with open(outputreferencefasta, "w") as refoutput:
            refoutput.write( str(">" + referencegene + "\n" ))
            refoutput.write(referenceseq)


        mafftoutput = os.path.join(OutputMAFFT_SNPs_Outputs, geneReformatted + "_" + patientNum + "_" + CC + "_" + Cluster1 + "_" + Cluster2 + "_" + Phenotype  +"_MAFFT.aln" )
        snpsoutput = os.path.join(OutputMAFFT_SNPs_Outputs, geneReformatted + "_" + patientNum + "_" + CC + "_" + Cluster1 + "_" + Cluster2 + "_" + Phenotype  +"_snpsites" )


        mafftcommand = "mafft --localpair --addfragments " + outputfastapath +" " + outputreferencefasta + " > " + mafftoutput + "\n"
        mafftsedcommand = "sed -i 's/\\" + "-/" + "\*/g' " + mafftoutput
        snpsitescommand = "snp-sites -v -o " + snpsoutput + " " + mafftoutput + "\n"

        comparisonsMAFFTcommands.append(mafftcommand)
        comparisonsSEDcommands.append(mafftsedcommand)
        comparisonsSNPcommands.append(snpsitescommand)
        # > /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_patient_159/Patient159/GeneFastas/MAFFT/Mafft_group_3858

        if numUniquegenes > 0 :
            mafftoutputfilename = "Patient_" + patientNum + "_" + CC + "_" + Cluster1 + "_" + Cluster2 + "_" + Phenotype + "_MafftCalls.sh"
            snpsitesoutputfilename = "Patient_" + patientNum + "_" + CC + "_" + Cluster1 + "_" + Cluster2 + "_" + Phenotype + "_SNPsites.sh"
            mafftSEDfilename = ("Patient_" + patientNum + "_" + CC + "_" + Cluster1 + "_" + Cluster2 + "_" + Phenotype + "_SEDcalls.sh")
            mafftoutputfilepath = os.path.join(OutputMAFFT_SNPsPath,str(mafftoutputfilename) )
            SNPs_outputfilepath = os.path.join(OutputMAFFT_SNPsPath,snpsitesoutputfilename )
            SEDfilespath = os.path.join(OutputMAFFT_SNPsPath,mafftSEDfilename )

            headerstring="#!/bin/bash\n"
            activatestring="mamba ~/mambaforge/bin/activate IntraPatientSNPsEnv \n"

            mafftcalls = open(mafftoutputfilepath, "w")
            mafftcalls.write(headerstring)
            mafftcalls.write(activatestring)

            snpcalls = open(SNPs_outputfilepath, "w")
            snpcalls.write(headerstring)
            snpcalls.write(activatestring)

            sedcalls = open(SEDfilespath, "w")
            sedcalls.write(headerstring)
            sedcalls.write(activatestring)

            for i in range(len(comparisonsMAFFTcommands)):
                mafftcalls.write(comparisonsMAFFTcommands[i])
                mafftcalls.write("\n")

                snpcalls.write(comparisonsSNPcommands[i])
                snpcalls.write("\n")

                sedcalls.write(comparisonsSEDcommands[i])
                sedcalls.write("\n")


        else:
            print("No mismatching core genes for this comparison")
