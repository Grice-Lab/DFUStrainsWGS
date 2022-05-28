# Amy Campbell
# 2022
# Plotting isolates from Patient 141

library("ape")
library("ggtree")
library("dplyr")
library("ggplot2")
library("ggtreeExtra")


TreeFilePath = "data/Patient141tree/RAxML_bestTree.RaxMLTree141Isolates"
UTDPhenotypes=read.csv("data/staphyloxanthin_paper_data.csv")
StaphIsolateDORNs$DORN=paste0("DORN", StaphIsolateDORNs$Doern.lab.bank.)

TreeObject = ggtree::read.tree(TreeFilePath)
RootedTree = ape::root(TreeObject, "CC1_MW2",  resolve.root=T)
ggtree(RootedTree)

# Adjust Outgroup length
TreePlot <- ggtree::ggtree(RootedTree, size=.5, layout = "rectangular") + geom_tiplab(size=2)

# Scale down the outgroup by 1/10
TreePlot$data[TreePlot$data$label %in% c("CC1_MW2"), "x"] = (TreePlot$data[TreePlot$data$label %in% c("CC1_MW2"), "x"])*.01

TreePlot + geom_treescale(label = "Estimated # substitutions per site (across 2442 core genes)") 

