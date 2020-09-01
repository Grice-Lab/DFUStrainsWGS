# Amy Campbell
# April 2020
# Makes a tab-delimited record, for each 'hypothetical protein group' defined by
# Prokka and Roary, of the 'Gene' (usually group_i), the header of first column
# it identifies which is not NA (AKA the first genome that contains the protein)
# the value of that column in that row (the locus tag), the same but for the
# second non NA column, third...etc.


# For a given row(gene) vector:
# Count non NAs in the

import pandas as pd

# For a given row, takes the genome-specific gene names for a given gene group
# from the first five genomes which contain it(if there are more than 5, that is.
# Otherwise, it takes as many as it has)
# Returns a little data frame which will go into the behemoth dataframe we'll use to
# map most likely genes back to 'hypothetical protein' annotations.
def return_rows(pandarow):
    General_Gene = (pandarow['Gene'])
    pandarow = pandarow.iloc[14: ]

    # Count # of non-NA columns representing genomes (so cols 14 and later)
    genomecount = int(pandarow.transpose().count())
    cols=pandarow.dropna().transpose()
    cols=cols.reset_index(drop=False)
    cols.columns=['Genome', 'GeneID']
    cols['CommonGeneName'] = General_Gene
    return(cols)
    #print(genomecount)
    # if genomecount <= 100:
    #     cols = pandarow.dropna().transpose()
    #     cols = cols.reset_index(drop=False)
    #     cols.columns = ['Genome', 'GeneID']
    #     cols['CommonGeneName'] = General_Gene
    #     return(cols)
    # else:
    #     cols = (pandarow.dropna())[:100]
    #     cols = cols.transpose()
    #     cols = cols.reset_index(drop=False)
    #     cols.columns = ['Genome', 'GeneID']
    #     cols['CommonGeneName'] = General_Gene
    #     return(cols)

# Read in the dataframe
BigDF = pd.read_csv("gene_presence_absence.csv", dtype=str)

# Filter the whole dataframe to rows where the 'Gene' column's value contains "group_"
HypDF = BigDF[BigDF['Gene'].str.contains("group_")]

# Apply the 'return_rows()' function to each row of the filtered DF.
obj = list(HypDF.apply(return_rows, axis=1))
pando = pd.concat(obj)
pando = pando.reset_index(drop=True)
pando.to_csv('Genes_to_Look_Up.csv')
