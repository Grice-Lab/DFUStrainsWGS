# Amy Campbell
# Python3 script

from Bio import SeqIO
import os
# testfasta = "/Users/amycampbell/Documents/PGAP/StaphGenomesAnnotation/DORN882/DORN882_pgap.fasta"
# newfasta = "/Users/amycampbell/Documents/PGAP/StaphGenomesAnnotation/DORN882/DORN882_TEST.fasta"

folderpath="/Users/amycampbell/Documents/PGAP/StaphGenomesAnnotation/"

pathlist = os.listdir(folderpath)

bigshellscript = open("/Users/amycampbell/Documents/PGAP/StaphGenomesHPC.sh", "w")

#print(pathlist)
for foldername in pathlist:
    fastafolder = os.path.join(folderpath, foldername)
    fastapath = os.path.join(fastafolder, (str(foldername)+"_pgap.fasta"))

    newrecords=[]

    #(1) Fix the contig header & save it anew
    ##########################################
    for record in SeqIO.parse(fastapath, "fasta"):
        otherstring = str(record.description)
        contigid = str(record.id)
        if str("circular=true") in otherstring:
            newID = "contig" +str(contigid) + " [topology=circular]"
        else:
            newID = "contig" + str(contigid)
        if len(record.seq) >= 200:
            newrecords.append(SeqIO.SeqRecord(record.seq, id=newID, description=""))
        else:
            print("Removed " + str(newID) + " length" + str(len(record.seq)))


    SeqIO.write(newrecords, fastapath, "fasta")
    #
    # # (2) Write the input.yml and the submol.yml for each one and add their calls to the shell script
    # #################################################################################################
    configfile1 = os.path.join(fastafolder, (str(foldername)+"_input.yml"))
    configfile2 = os.path.join(fastafolder, (str(foldername)+"_submol.yml"))

    openconfig1 = open(configfile1, "w")

    openconfig1.write("fasta:\n")
    openconfig1.write(" class: File\n")
    openconfig1.write(" location: " + str(foldername)+"_pgap.fasta\n")
    openconfig1.write("submol:\n")

    openconfig1.write(" class: File\n")
    openconfig1.write(" location: " +(str(foldername)+"_submol.yml\n"))

    openconfig1.close()


    openconfig2 = open(configfile2, "w")
    openconfig2.write("organism:\n")
    openconfig2.write("    genus_species: Staphylococcus aureus")
    openconfig2.close()




    outputstring = "./pgap.py -r -o staph_output/"+ str(foldername)+ "output StaphGenomesAnnotation/" + foldername + "/" +str(foldername)+"_input.yml"
    bigshellscript.write(outputstring)
    bigshellscript.write('\n')
    bigshellscript.write('\n')

bigshellscript.close()
