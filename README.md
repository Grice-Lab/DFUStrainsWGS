# DFU100 *S. aureus* Isolates
Full list of commands for the 2020 Assembly of 221 S. aureus isolates from the DFU100 cohort (after removal of non-*S. aureus* isolates), as sequenced on the Illumina HiSeq in late 2019, and the subsequent coverage & gene content analyses on these assembled genomes / selection of a subset of these isolates to perform Oxford nanopore long read sequencing on. Some 'assembly' portions of this README are shared with that of [EAGenomeAssembly](https://github.com/Grice-Lab/EAGenomeAssembly). Steps 1-4 were performed on the 223 (presumed\*) *S. aureus* isolates at the same time as they were performed on some other submitted isolates. 

\* One of the 223 DFU genomes turned out to be *S. xylosis*, and another turned out to be *S. simulans* 
## Conda Environments
Conda environment configuration files for this project are located [here](https://github.com/Grice-Lab/DFUStrainsWGS/tree/master/environments)
- **EAGenomeEnv**
- **BlastEnv**
- **BowtieEnv**
- **QuastEnv**
- **R_envir**
- **TreeEnv**
- **TrimmingEnvironment**
- **updated_annotation_env**

## **Assembly and Quality Control**
### 1. ) FastQC on untrimmed, demultiplexed Illumina reads 
Here, we run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the raw, demultiplexed reads to assess the quality of the reads prior to trimming adapters or trimming for length. ***Outputs quality report***


#### Commands
Shell script runs FastQC in paired end mode with default settings from **EAGenomeEnv**.
>  [RunFastQC_PreTrim.sh](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/FastQC/RunFastQC_PreTrim.sh)

### 2. ) Trim ends with TrimGalore
Here, we run [TrimGalore](https://github.com/FelixKrueger/TrimGalore) to remove adapter sequences from the raw illumina reads and trim noisy ends observed in the fastqc step. **Make_TrimGalore_shells.py** builds shell scripts which contain commands to run TrimGalore  **TrimmingEnvironment** in paired-end mode, set to stringency of 10, with length cutoff of 70bp, and set to clip 10bp off of the 5' end of each read. These shells also make use of the --nextera flag to specifically remove the nextera transposase sequence. 
***Outputs trimmed reads***
#### Commands

From TrimmingEnvironment , make shell scripts.

> OutPutShells="/home/acampbe//EAGenomeAssembly/TrimGaloreShells"  
> OutPutTrimGalore="/project/grice/storage/HiSeq/WGS/HiSeq_19/TrimmedFastqs_TrimGalore"  
> CondaPath="/home/acampbe/software/miniconda3/bin/activate"  
> 
> python3 Make_TrimGalore_shells.py --inputdirname $InputDirectory --outputdirshells $OutPutShells --outputdir_trimgalore $OutPutTrimGalore --conda_activatepath $CondaPath --nshells 10

The following [shell scripts](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/TrimGaloreShells) run TrimGalore on each isolate:

> Run_TrimGalore_0.sh
> Run_TrimGalore_1.sh  
> Run_TrimGalore_2.sh  
> Run_TrimGalore_3.sh  
> Run_TrimGalore_4.sh  
> Run_TrimGalore_5.sh  
> Run_TrimGalore_6.sh  
> Run_TrimGalore_7.sh  
> Run_TrimGalore_8.sh  
> Run_TrimGalore_9.sh

### 3. FastQC on trimmed Illumina reads 
After running the trimming shells, show that TrimGalore actually *improved* the quality of our reads for input into the assembly  step by running FastQC again. ***Outputs updated quality report***
#### Commands
> [RunFastQC_PostTrim.sh](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/FastQC/RunFastQC_PostTrim.sh)

### 4. SPAdes assembly with Unicycler 
Run [SPAdes](http://spades.bioinf.spbau.ru/release3.11.1/manual.html) using [Unicycler](https://github.com/rrwick/Unicycler)'s wrapper for it on the illumina short reads only. ***Outputs assembled contigs***

#### Commands

 A python script generates [20 shellscripts](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/UnicyclerIlluminaShells) which, from EAGenomeEnv, calls unicycler on the paired-end short reads with default settings (running SPAdes on a range of k-mers)
> python3 [Make_UnicyclerIlluminaShells.py](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/MakeShells/Make_UnicyclerIlluminaShells.py) --inputdirname /project/grice/storage/HiSeq/WGS/HiSeq_19/TrimmedFastqs_TrimGalore --outputdirshells UnicyclerIlluminaShells --outputdir_unicycler /project/grice/storage/HiSeq/WGS/HiSeq_19/UnicyclerAssemble --conda_activatepath /home/acampbe/software/miniconda3/bin/activate --nshells 20

For 1...20, run 
> Run_UnicyclerSpades_<_>.sh


### 5. Align trimmed reads directly to the SPAdes assembly for each isolate to evaluate coverage and depth. 
Using **BowtieEnv**, the following shell scripts perform a [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) alignment of the trimmed reads to the SPAdes-generated assemblies from step 4, then use samtools to extract depth-by-contig information from these alignments.  
***Outputs bowtie alignments of trimmed reads against each isolate, depth, length for each contig in each isolate (DFU_Covg_stats.tsv)***
#### Commands  
Run Bowtie2, Samtools, and Galaxy tools
> [Run_Bowtie_DORN.sh](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/CoverageAnalyses/DFU_Coverage/Run_Bowtie_DORN.sh)  
> [Coverage_Stats_Contigs_DFU.sh](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/CoverageAnalyses/DFU_Coverage/Coverage_Stats_Contigs_DFU.sh)

Make output .csv file of depths and lengths by contig of each isolate
> [Merge_Coverage_Stats_DFU.sh](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/CoverageAnalyses/DFU_Coverage/Merge_Coverage_Stats_DFU.sh)


### 6. Filter contigs by length and depth
After the initial contig assembly, in order to make sure we're making conclusions about gene content and SNP-level relationships between isolates based on a real representation of each genome and not based on artifacts of either contamination or misassembly. While the latter is more difficult to detect in the absence of long reads against which we can scaffold, we should at least filter out contigs which are too short, of too low coverage, or map to another organism than *S. aureus* by BLAST. ***Outputs record of retained versus removed contigs, cleaned contigs for core genome analysis and tree building***

We first sort the contigs into three categories:   
* **Remove**: Contigs which are < 500bp in length or <10X depth
* **Follow up**: Contigs which have a depth less outside of the **(Median Contig Depth for the Isolate)** ± **1.5 x (Depth IQR)** range
* **Keep**: Contigs which fit none of the above criteria 

Then, we search the contigs labeled for "follow up" against the BLAST nt database. We note any follow-up contig (high- or low-coverage outlier) for which the best hit (by bitscore) is not of the *Staphylococcus aureus* species (this turned out to only be one, low-coverage outlier contig, which was circular and mapped closely to known *S. aureus* plasmids, and we won't remove it). 


#### Commands

Run R script to classify contigs' keep vs. followup vs. drop status based on depth and length
> [Contig_Depth_Analyses.Rmd](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/CoverageAnalyses/DFU_Coverage/Contig_Depth_Analyses.Rmd)
> 
Run a shell calling the python script to write .fastas based on the depth-based contig classification
> [SortContigsDepth.sh](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/CoverageAnalyses/DFU_Coverage/SortContigsDepth.sh)

Blastn query for every contig classified for 'followup,'

> [BlastFollowups.sh](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/CoverageAnalyses/DFU_Coverage/BlastFollowups.sh)

Add back the "follow-up" contigs that don't need to be removed (none of them did, in fact)
This shell script uses **Add_SA_contigs_Back.R** to parse blast output, make fasta of 'added back' contigs 
> [ReadBlastOutput.sh](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/CoverageAnalyses/DFU_Coverage/ReadBlastOutput.sh)  
> [CombineFilteredContigs.sh](https://github.com/Grice-Lab/EAGenomeAssembly/tree/master/CoverageAnalyses/DFU_Coverage/CombineFilteredContigs.sh)

## **Phylogeny**
From the cleaned assemblies, we build a SNP-based, maximum likelihood phylogeny of the 221 *S. aureus* isolates with the following goals:
- Visualize and interpret the evolutionary relationships between the different isolates observed in the DFU100 cohort.
- Identify subject- and healing outcome-specific clades of *S. aureus*
- Interpret the isolate genomes in the context of known clinical strains of *S. aureus* by estimating the positions of each reference genome *a posteriori*

### 1. Automated prokka annotation of each genome
Shellscripts contained [here](https://github.com/Grice-Lab/DFUStrainsWGS/tree/master/Phylogeny/ProkkaShells) call Prokka on default settings to predict gene coordinates and protein sequence. ***outputs .gff annotation files, gene coordinate maps, and .faa protein translations for each genome***  
#### Commands
For 1...9 run:
> Run_Prokka<>.sh  

Run Prokka on the selected *S. aureus* reference genomes and outgroup *S. epidermidis*
> Run_Prokka_references.sh  
### 2. Run [Roary](https://github.com/sanger-pathogens/Roary) 
Define core & accessory genomes for the set of isolates and build multiple alignment between core genomes. 
The [pplacer](https://github.com/matsen/pplacer) algorithm, which will place the reference genomes on the maximum likelihood phylogeny of isolates, requires an alignment between all isolates and the references as input. Therefore, we will also run roary (and its wrapper for [PRANK](http://wasabiapp.org/software/prank/)) on the isolates in aggregate with the reference genomes. ***outputs core and accessory genomes for each isolate, probabilistic multiple alignment of core genomes for the set with and without the references/outgroup included, SNP distance matrix for each alignment***

#### Commands

>[Run_Roary.sh](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/Run_Roary.sh)  
>[Run_Roary_References.sh](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/Run_Roary_References.sh)  
>[Run_SNPdists.sh](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/Run_SNPdists.sh)
 ### 3. Build maximum likelihood tree 
Run [RaxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) to build a maximum likelihood, SNP-based tree from the core genome alignment of the 221 isolates with a random seed set, using the GTR-gamma nucleotide subsitution model. ***outputs a maximum likelihood tree in \.newick format, as well as a text file containing the parameters of the final ML model.***
#### Commands
> [Run_RaxMLCore.sh](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/Run_RaxMLCore.sh)

### 4. Place reference *S. aureus* genomes and *S. epidermidis* outgroup onto the ML tree *a posteriori*  
Pplacer estimates the maximum posterior probable location of each reference genome on the maximum likelihood tree  output by RaxML. **Note:** A silly versioning issue(pplacer was developed on a different version of RaxML than I used) means that I had to manually remove the "Partition: 0 with name: No Name Provided” string from the RAxML parameter file(input with the '-s' flag into Pplacer). ***outputs a .jplace file containing the mappings of each reference genome to the tree, and a .newick representation of this mapping thanks to [tog](https://matsen.github.io/pplacer/generated_rst/guppy_tog.html)***  

#### Commands
> [RunPplacer.sh](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/RunPplacer.sh)

### 5. Visualize RaxML/Pplacer Tree
This R script visualizes the ML tree of isolates with references placed on it, collapsing clades of closely related genomes, recording the descendants contained in each collapsed clade. ***outputs plots for tree visualization as well as a mapping of the closest reference genome by SNP distance to each of the isolates' core genomes***
#### Commands
> Rscript [TreeMapping_ShortReads_CleanedContigs.R](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/TreeMapping_ShortReads_CleanedContigs.R)
## **Final assembly quality and genome characteristics**

- Get sequencing depth and breadth of coverage estimates for each isolate
- Summarize confirmed circular contigs (e.g. plasmids) present in each isolate's genome
- Identify regions of each genome which differ from reference genomes most closely related by core SNP distance 
- Select a subset of the 221 *S. aureus* isolates to perform Oxford Nanopore long read sequencing on. 


### 1. Calculate average depth and breadth of coverage for each assembly

Here, we calculate ***depth*** as the average # of bases mapped to each position on an assembly, which we extract from Bowtie2 alignments of each isolate's trimmed reads to its assembly. We calculate ***breadth of coverage*** as the proportion of positions along a reference genome that each isolate's reads map to with at least 10X depth. We extract this information from Bowtie2 alignments of each isolate's trimmed reads to the closest (by core SNP distance). The "closest" reference on the tree to each isolate was output in step 2 of the previous section. 

### Commands
This script runs the bowtie alignments for both cases (isolate's reads against its assembly, isolate's reads against its closest reference genome) ***outputs one text file with the average depth per isolate, another text file with the base coverage count for each isolate to its reference, as well as reference length***
> [Coverage_Stats_CleanedContigs_DFU.sh](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/CoverageAnalyses/DFU_Coverage/Coverage_Stats_CleanedContigs_DFU.sh)

### 1. Calculate assembly quality statistics for each assembly
Calculate statistics such as total # contigs, N50, L50 on the cleaned assemblies  using [QUAST](http://bioinf.spbau.ru/quast). 

> [Run_Quast_DORN.sh](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/CoverageAnalyses/DFU_Coverage/Run_Quast_DORN.sh)

### 2. Count the # of circularized plasmids identified in each assembly.

> [Count_Plasmids.sh](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/CoverageAnalyses/DFU_Coverage/Count_Plasmids.sh)

### 3. Identify regions of each isolate's genome that don't map to the nearest reference  
Run MiniMap2 pairwise alignments between each isolate's assembly and the reference genome on the tree it's closest to (by SNP distance). This enables us to identify areas of an isolate that aren't represented in its closest reference genome.  ***Outputs tabular summary of pairwise alignment***

### Commands
> [Run_MiniMap2.sh](https://github.com/Grice-Lab/EAGenomeAssembly/blob/master/CoverageAnalyses/DFU_Coverage/Run_MiniMap2.sh)
### 4. Aggregate statistics from steps 1-3 to select representative isolates from each cluster in step 1 for ONP sequencing

> Rscript [AggregateGenomeStats.R](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/AggregateGenomeStats.R)

### 5. Identify "clusters" of isolate genomes which share >99.99% core genome identity. Select representative isolates for each cluster based on the the output of Step 5. 
The **SubsettingCleanedAssemblies_forONP.R** script constructs groups of isolates from the same subject that all differ by less than 118 SNPs in their core genomes (.01% the length of the core genome alignment between all isolates). It selects 87 representative isolates(to fill a 96-well plate) from the set based on the following criteria: 
- All subjects represented in the set of 221 isolates should have at least one isolate sequenced with ONP
- Any isolate which differed by more than >118 SNPs from all other isolates should be sequenced with ONP (44 isolates total)
- At least one isolate should be sequenced from each core genome-based cluster. The representative should be chosen with the following "prioritized" criteria:
  * High enough depth so we have at least one really good hybrid assembly from that ONP scaffold
  * Imperfect coverage to the closest reference genome (step 1)
  * Prioritize isolates with excess unmapped sequence to their closest reference by pairwise alignment (step 3)
  * Fewer plasmids which have already been circularized using short read sequencing/ plasmidSPAdes
  * Fragmented assemblies with many shorter contigs  
- While representatives for the 11 largest clusters were chosen manually (>7 isolates each) based on the above criteria and plots of breadth vs. coverage, representatives for the 12th:42nd clusters were chosen in an automated fashion with the following nested criteria:  
    * Maximum total unmapped sequence to closest reference genome. 
    * Most fragmented assembly, based on minimum N50 score, if unmapped lengths equal
    * Minimum number of circularized plasmids if N50 equal
 

>[SubsettingCleanedAssemblies_forONP.R](https://github.com/Grice-Lab/DFUStrainsWGS/blob/master/Phylogeny/SubsettingCleanedAssemblies_forONP.R)







 


