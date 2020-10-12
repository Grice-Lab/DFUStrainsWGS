## 1.)  See whether the entire cassette is on the same contig
**INPUT:**
 - *--Genomes_Exclude* List of genomes to exclude(DORN1340,DORN1473,DORN672,DORN691,DORN701 by default since they don't have complete versions of crtOPQMN genes)
 - *--geneNames* : "group_2428,crtP,crtQ,crtM,crtN" by default; these are the shared ID's spit out by **Roary** for the crt cassette genes
 - *--gffDirPath* : path to a folder containing all of the .gff files output by Prokka and used for the run of Roary
 - *--PresenceAbsencePath*: path to the Roary-generated 'gene_presence_absence.csv' file

**OUTPUT**:
 - *XanthinContigInfo.csv*, which indicates, for each genome considered, whether the genes are all on the same contig, which contig they're on, and the beginning/ending loci for the cassette in the genome. 

**COMMANDS**:
> conda activate xanthinEnv
> python3 Summarize_crtOPQMN.py

## 2.) Identify contig & loci containing promoter sequence for the operon 
**INPUT:**
 - *xanthinpromoter.fasta* : sequence of the promoter from  Staphylococcus aureus subsp. aureus strain Dresden-275757 chromosome (CP054876.1) as identified by a BLASTN search for the primers used in Pelz et al.'s 2005 [study of the staphyloxanthin operon structure](http://doi.org/10.1074/jbc.M505070200) . 
 - Cleaned contigs of the genome assemblies (folder containing assemblies in .fasta form)

**OUTPUT**: 
- a folder called *XanthinPromoterBlastResults/* containing .tab output for the blast searches of the promoter reference sequence against each genome assembly. 
- These are aggregated into *AllPromoterSequences.tab*, which is a .tsv containing the fields: GenomeID, queryID, contigID within the Genome, PctIdentity, Length of Hit, Starting position in contig, ending position, e-value, bit score. 

**COMMANDS**:
> bsub -e promotersearch.e -o promotersearch.o sh BlastPromoter.sh
> bsub -e parse.e -o parse.o sh ParseBlastResults.sh

