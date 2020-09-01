# Amy Campbell
# 08/2020

# Parses the messy bash output from two files made with SAMTOOLS

# breadthcounts.txt
# example entry:
# DORN1000
# 2750256 --> # bases > 10X depth
# 106.945 294251724 --> IGNORE THIS LINE OK
# 2820462 --> reference length
# So, breadth of covg is 2750256/2820462 or .975



# 2750256 = # bases with >10x depth when you align the isolate's trimmed reads
# to the reference genome
#  =
Breadthoutput = open("BreadthData.txt", "w")
Breadthoutput.write("DORN\tBreadthOfCovg\n")

breadths = open("/Users/amycampbell/Desktop/Club_Grice/Club_Grice/EAGenomeAssembly/CoverageAnalyses/DORN_CoverageAnalyses/breadthcounts.txt","r")
breadthlist = breadths.readlines()
breadths.close()

for b in range(0, len(breadthlist), 4):

    DORN= str(breadthlist[b])
    DORN= DORN.strip("\n")
    total10X = int(breadthlist[(b+1)])
    totallength = int(breadthlist[(b+3)])
    breadthofcovg= total10X/totallength
    print(DORN)
    print(breadthofcovg)
    Breadthoutput.write( str(DORN) + "\t" + str(breadthofcovg) + "\n")
Breadthoutput.close()





DepthOutput = open("DepthData.txt", "w")
depths = open("/Users/amycampbell/Desktop/Club_Grice/Club_Grice/EAGenomeAssembly/CoverageAnalyses/DORN_CoverageAnalyses/DORNDepths.txt", "r")

DepthOutput.write("DORN\tDepth\n")

depthlist = depths.readlines()
depths.close()

for d in range(0, len(depthlist), 3):
    DORN= str(depthlist[d])
    DORN= DORN.strip(".sam\n")

    depthsline = depthlist[(d+2)]
    depthslinelist = depthsline.split()
    averagedepth = str(depthslinelist[0])
    DepthOutput.write(str(DORN + "\t" + averagedepth + "\n"))


DepthOutput.close()
