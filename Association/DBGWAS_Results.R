# Amy Campbell
# December 2020 
# Mapping DBGWAS Results onto the tree

# Slogging through results from DBGWAS
# filtering to components containing unitigs significant at q value < .05 level 

##################################
# READ IN PHYLOGENY FOR REFERENCE
##################################
setwd('/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/')
library("treeio")
library("ggtree")
library("dplyr")
library("ggplot2")
library("dplyr")
library("ape")
library("tibble")

SubjTime = function(row){
  if(is.na(row["DORN"])){
    return(NA)
  }else if(is.na(row["patient_id"])){
    return(row["DORN"])
  }else{
    return(paste(row["patient_id"], stringr::str_replace(row["visit"], " ", ""), sep="_"))
  }}

Isolates = "data/DFU_Staph_aureus_isolates.csv"
CoreTree_PosteriorReferences = "data/core_gene_alignment.newick"
ViewTree=treeio::read.jplace("data/core_gene_alignment.jplace")

IsolatesInfo = read.csv(Isolates)
IsolatesInfo$DORN = paste0("DORN", IsolatesInfo$Doern.lab.bank.)
IsolatesInfo$subject_timepoint = paste(IsolatesInfo$patient_id, IsolatesInfo$visit,sep="_")
IsolatesInfo$Genome2=paste0("DORN", IsolatesInfo$Doern.lab.bank.)
IsolatesInfo$Genome1=paste0("DORN", IsolatesInfo$Doern.lab.bank.)
ReferenceNames=c("SA_502A","S_epidermidis",
                 "USA100_AR465", "USA300_FPR3757",
                 "USA400_051", "CC8_NCTC8325",
                 "CC8_Newman", "CC1_MSSA476",
                 "CC1_MW2", "CC22_HO",
                 "CC30_MRSA252", "CC398",
                 "CC398_ATCC6538","CC5_Mu50",
                 "CC5_N315","CC72_CN1", "SA_CFSAN007883", "SA_MCRF184", "SA_AR464", "SA_UP_1150")
PPlacerTree = ggtree::read.tree(CoreTree_PosteriorReferences)
apeRooted = ape::root(PPlacerTree,"S_epidermidis", resolve.root=T)

# Adjust Outgroup length
apeplot <- ggtree::ggtree(apeRooted, size=.25, layout = "rectangular") + geom_tiplab(size=2)
apeplotbackup = apeplot
apeplot$data$DORN = apeplot$data$label
apeplot$data[apeplot$data$label %in% c("S_epidermidis"), "x"] = mean(apeplot$data$x)
apeplot$data$DORN = if_else(apeplot$data$DORN == "CC398_ATCC6538", "ATCC6538", apeplot$data$DORN )
apeplot$data$label = apeplot$data$DORN

############################
# STAPHYLOXANTHIN PHENOTYPE
############################
Unitigs_By_Isolate_Path = "/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/Results/DBGWAS_Output_WithAnnotations/step1/bugwas_input.all_rows.binary"
Results = "/Users/amycampbell/Box/GRICE\ LAB\ SHARE/Current\ lab\ members/Amy/StaphStrainAnalysis/Results/DBGWAS_Output_WithAnnotations/textualOutput/all_comps_nodes_info.tsv"

UnitigDist = read.table(Unitigs_By_Isolate_Path,header = T )
UnitigDist$node = paste0("n", UnitigDist$ps)


# Comp_33 5 Significant (FDR-adjusted q value <.05) unitigs mapped to FCOIFMIJ_01968_aroK
################################################################################
aroK_nodes = c("n118121", "n118247", "n291792", "n233670", "n72782")
aroK_dist = UnitigDist %>% filter (node %in% aroK_nodes)
aroK_dist_ref = t(aroK_dist[,2:(ncol(aroK_dist) - 1)])
colnames(aroK_dist_ref) = aroK_dist$node
View(aroK_dist_ref)

aroK_dist_ref$sum_variants = aroK_dist_ref$n72782 + aroK_dist_ref$n118121 + aroK_dist_ref$n118247 + aroK_dist_ref$n233670 + aroK_dist_ref$n291792
aroK_dist_ref = data.frame(aroK_dist_ref)
aroK_dist_ref$DORN = rownames(aroK_dist_ref)
aroKplot = apeplot
aroKplot$data = aroKplot$data %>% left_join(aroK_dist_ref, by="DORN")
aroKtree = ggtree(aroKplot$data) + geom_point(aes(color=(n233670>0),shape=(n233670>0), size=10))+ scale_color_manual(values=c("NA","darkgoldenrod2")) + geom_tiplab()
ggsave(aroKtree, file="Node233670Tree.pdf",width=40,height=48)
# Tree where tips are colored by # of non-mutually exclusive K-mer variants a variant contains 
aroKtree_any = ggtree(aroKplot$data) + geom_point(aes(color=as.factor(sum_variants),shape=(sum_variants>0), size=10))+ scale_color_manual(values=c("NA", RColorBrewer::brewer.pal(name="YlOrRd", n=5))) + geom_tiplab()
ggsave(aroKtree, file="aroK_Variant_Tree.pdf",width=40,height=48)


ggtree(aroKplot$data) + geom_point(aes(color=as.factor(sum_variants),shape=(sum_variants>0), size=10))+ scale_color_manual(values=c("NA", RColorBrewer::brewer.pal(name="YlOrRd", n=5))) + geom_tiplab()
ggtree(apeplot)


# 3 Significant ((FDR-adjusted q value <.05) unitigs ADJACENT to those
# mapped to FCOIFMIJ_00620_group_8709 (an acyltransferase family protein; core gene in this set)
#######################################################################
group_8709_adjacent_nodes = c("n235801", "n33647", "n790")
group_8709_dist = UnitigDist %>% filter (node %in% group_8709_adjacent_nodes)
group_8709_dist_ref = data.frame(t(group_8709_dist[, 2:(ncol(group_8709_dist)-1)]))
colnames(group_8709_dist_ref)= group_8709_dist$node
group_8709_dist_ref$DORN = rownames(group_8709_dist_ref)
group_8709_dist_ref$sum_variants = group_8709_dist_ref$n790 + group_8709_dist_ref$n33647 + group_8709_dist_ref$n235801

group_8709plot = apeplot
group_8709plot$data = group_8709plot$data %>% left_join(group_8709_dist_ref, by="DORN")
group_8709tree = ggtree(group_8709plot$data) + geom_point(aes(color=as.factor(sum_variants),shape=(sum_variants>0), size=10))+ scale_color_manual(values=c("NA", RColorBrewer::brewer.pal(name="YlOrRd", n=3))) + geom_tiplab()
ggsave(group_8709tree,width=40,height=45, file="group8709" )

# Comp_45 One unitig mapped to and one adjacent to isdG_1 (group_1221, Heme oxygenase (staphylobilin-producing))
#############################################################################################
isdg_1_adjacent = c("n215989",
                    "n289588")

# Comp_90 4 unitigs adjacent to yagU https://bmcmicrobiol-biomedcentral-com.proxy.library.upenn.edu/articles/10.1186/1471-2180-6-89
#############################################################################################
yagUadjacent = c("n69288",
  "n217233",
  "n286547",
  "n69289")

# Comp_181 2 unitigs adjacent to hdfR_1
#############################################################################################
c("n261342",
"n269461")

# # Comp_24 3 Significant unitigs mapped to regions between those annotated to group_8709
# and those mapped to slyA_1
###############################################################################################
c("n790,
n235801,
n33647")

# Comp_206 4 unitigs with low minor allele frequency annotated to group_2585 or GTP cyclohydrolase 1 type 2
###################################################################################################
c("n172254",
"n277465",
"n279704",
"n213212")


# Comp_135 8 unitigs (3 annotated to group_3098 /FCOIFMIJ_0699; 5 annotated to isdE
###################################################################################################
c("n235355",
  "n45467",
  "n45468",
  "n211744",
  "n291639",
  "n240183",
  "n245639",
  "n61643")


# Comp_30 16 significant unitig
###################################################################################################
# n144358, n279428, n243027, n13166 close to a region annotated to group_4701, group_4883, group_5062, group_4942
# n117965 between regions annotated to the above and regions annotated to group_4702
# n103950, n185153, n309678, n55873, n243409 annotated to group_4702 
# n153107 and n192918 annotated to group_542, group_543, group_544, group_545
# n7041 annotated to group_4703, group_4704, group_545(larger e-value)
# n57034 not annotated but close to 7041
# n109758 annotated to group_3569, dut, group_5464, group_6000
# n99239 annotated to group_169, group_1141, group_165

comp30_nodes = c("n7041",
  "n13166",
  "n55873",
  "n57034",
  "n99239",
  "n103950",
  "n109758",
  "n117965",
  'n144358',
  "n153107",
  "n185153",
  "n192918",
  "n243027",
  "n243409",
  "n279428",
  "n309678")




