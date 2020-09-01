# Amy Campbell 02/2020
# Assign AGR types to genomes based on output of the tblastN search of AGR-D AIP sequences against the whole nucleotide sequences of each assembly 
# To be used by ARM for quorum phenotyping assays

setwd("/Users/amycampbell/Desktop/Club_Grice/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts")

library("dplyr")
library("ggplot2")

# Blast results are from TBlastN for the AIP (AGR-D) sequences contained in literature-assigned examples of each AGR type:

# Protein 'Queries' for each blast search 
#########################################
# TypeI: Strain V_2200, MNTLFNLFFDFITGILKNIGNIAAYSTCDFIMDEVEVPKELTQLHE
# TypeII: Strain 502A, MNTLVNMFFDFIIKLAKAIGIVGGVNACSSLFDEPKVPAELTNLYDK
# TypeIII: Strain MW2, MKKLLNKVIELLVDFFNSIGYRAAYINCDFLLDEAEVPKELTQLHE
# TypeIV: Strain NRS153, MNTLLNIFFDFITGVLKNIGNVASYSTCYFIMDEVEIPKELTQLHE


# Controls: NT sequences for the following
##########################################
# TypeI: GCF_000463055.1 (CN1)
# TypeII: GCF_001548415.1 (Strain MI)
# TypeIII: GCF_000772025.1 (strain FORC_001)
# TypeIV: GCF_000462955.1 (Strain 6850)
# AGR-KO: Took the MW2 genome, pseudorandomly mutated region around and including AGR-D's AIP (240 nts) making sure not to lose or introduce any stop codons 

# Combined_TblastN.tab just contains the concatenated results of the search described above
blastresults = read.csv("Combined_TblastN.tab", header=F, sep="\t")
colnames(blastresults) = c('Query', 'Subject', 'PctIdentity', 'PctQueryCoveredSubjectSequence', 'evalue','bitscore')


blastresults$Subject = sapply(blastresults$Subject, function(x) (strsplit(toString(x),'_'))[[1]][1])

# Look at E values for positive controls -- make sure 'correct' evalues are not neck and neck with next lowest e values

# Lowest e-value with type I, next lowest is type IV 3 orders of magnitude higher
blastresultsI = blastresults %>% filter(Subject=="TypeI")

# Lowest e-value with type II, next lowest is with type I, 7 orders of magnitude higher
blastresultsII = blastresults %>% filter(Subject=="TypeII")

# Lowest E value for type III, next lowest is 15 orders of magnitude higher 
blastresultsIII = blastresults %>% filter(Subject=="TypeIII")

# Lowest E value for type IV, next lowest is with type I at 4 orders of magnitude higher
blastresultsIV = blastresults %>% filter(Subject=="TypeIV")

# Get predicted type for each-- highest pct identity & lowest e value
#####################################################################
blastresults_Reduced = blastresults %>% group_by(Subject) %>% filter(evalue==min(evalue))

# Look at distribution of e values
ggplot(blastresults_Reduced, aes(x=log10(evalue))) + geom_histogram(binwidth=1) + ggtitle("E-values for AGR-D TBlastN Searches")


# Anything with an e value above, say, 10^-6should be listed as "unassigned"
# otherwise:
# if Query=="AGR_D_I_V2200", assign type I
# else if Query=="AGR_D_II_502A", assign type II 
# else if Query=="AGR_D_III_MW2", assign type III
# else if Query=="AGR_D_IV_NRS153", assign type IV
blastresults_Reduced = blastresults_Reduced %>% mutate(AGR_Type=case_when(
  (evalue > 1e-6)~ "Unassigned", 
  (evalue < 1e-6 && Query=="AGR_D_I_V2200") ~ "TypeI",
  (evalue < 1e-6 && Query=="AGR_D_II_502A") ~ "TypeII",
  (evalue < 1e-6 && Query=="AGR_D_III_MW2") ~ "TypeIII",
  (evalue < 1e-6 && Query=="AGR_D_IV_NRS153") ~"TypeIV"
))

write.csv( blastresults_Reduced,"blastResults_assignments.csv")
blastresults_Reduced_TypeI = blastresults_Reduced %>% filter(AGR_Type=="TypeI")
blastresults_Reduced_TypeI$Subject

blastresults_all_typeI = blastresults %>% filter(Subject %in% blastresults_Reduced_TypeI$Subject)
View(blastresults_all_typeI)
Sort_typeI_subjEvalue = blastresults_all_typeI %>% arrange(Subject, evalue)
write.csv(Sort_typeI_subjEvalue, "Evalues_TypeI_Assignments.csv")

blastresults_Reduced_TypeII = blastresults_Reduced %>% filter(AGR_Type=="TypeII")
blastresults_all_typeII = blastresults %>% filter(Subject %in% blastresults_Reduced_TypeII$Subject)
Sort_typeII_subjEvalue = blastresults_all_typeII %>% arrange(Subject, evalue)
VIEW(Sort_typeII_subjEvalu
    )

View(Sort_typeII_subjEvalue)

blastresults_Reduced_TypeIII = blastresults_Reduced %>% filter(AGR_Type=="TypeIII")
blastresults_all_typeIII = blastresults %>% filter(Subject %in% blastresults_Reduced_TypeIII$Subject)
Sort_typeIII_subjEvalue = blastresults_all_typeIII %>% arrange(Subject, evalue)

type4 = blastresults %>% filter(Subject=="TypeIV")

View(type4)
