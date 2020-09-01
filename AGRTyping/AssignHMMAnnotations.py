################################################################################
# Amy Campbell
# April 2020
# Reads the CSV of gene cluster ID's without annotations ('hypothetical')
# Iterates through each unique genome in the dataframe,
# sort its DF by GeneID's alphabetically
# Opens the <genome>_HMMerHits.tab
# Goes through the now alphabetical GeneIDs.
################################################################################

import pandas as pd
import os
# Set my local paths to the table produced by HypotheticalProteinGroupMap.py,
# Directory containing HMMer tables in <genomename>_HMMerHits.tab format

Lookup_Genes_Path = "/Users/amycampbell/Desktop/Club_Grice/Club_Grice/scripts/acampbe/DFU/scripts/isolates_analysis_scripts/Genes_to_Look_Up.csv"
HMMERtable_Path = "/Users/amycampbell/Desktop//HMMerTables"


def AssignAnnotation(geneID, genome):
    subset = genome.loc[genome['query']==geneID]

    if subset.empty:
        return ["None", 1.0000000000000]
    else:
        return (subset['target']).values[0], (subset['evalue'].values[0])


# print((pando['Genome'].unique()).shape)
LookUpGenes = pd.read_csv(Lookup_Genes_Path)
genome_String_list = (LookUpGenes['Genome'].unique())

NewDF = pd.DataFrame(columns=['target', 'accession', 'GeneID', 'accession_query', 'evalue', 'bitscore','genome','CommonGeneName'])
for genome_String in genome_String_list:
    genomelookup = LookUpGenes.loc[LookUpGenes['Genome']==genome_String]
    lookup_list =  genomelookup['GeneID']
    lookup_DF = genomelookup[['GeneID', 'CommonGeneName']]
    genomeDF = pd.read_csv(os.path.join(HMMERtable_Path, str( genome_String + "_HMMerHits.tab")),delim_whitespace=True, usecols=[0,1,2,3,4,5], skiprows=3, skipfooter=10, header=None, index_col=False)
    genomeDF.columns = ['target', 'accession', 'GeneID', 'accession_query', 'evalue', 'bitscore']
    genomeDF['genome'] = genome_String
    genomeDF = genomeDF.loc[genomeDF['GeneID'].isin(lookup_list.values)]
    indices = genomeDF.groupby(['GeneID'])['evalue'].transform(min) == genomeDF['evalue']
    genomeDF = genomeDF[indices]
    FinalDF = pd.merge(genomeDF, lookup_DF, on='GeneID')
    NewDF = pd.concat([NewDF,FinalDF])


    #for geneid in genomelookup['GeneID']:
    #    NewDF.at[Row_Indexer, 'HMMAnnotation'], NewDF.at[Row_Indexer, 'evalue'] =  AssignAnnotation(geneid, genomeDF)
    #    Row_Indexer = Row_Indexer + 1



NewDF.to_csv("FullMappingFile_HMMs.csv")
# While i < numRows_Genome:
#   query=genome[row i gene name]
#   iterate through rows in genome_HMMerHits tab until you find the query
#
#For the first Gene ID, find the first occurence of it.
