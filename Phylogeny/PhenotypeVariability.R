library(dplyr)
library(ggplot2)
setwd("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS")
phenotypes = read.csv("data/phenotype_variation_11.06.20.csv")

MinMaxNormalization = function(phenotype_col){
  minimum =  min(phenotype_col, na.rm=T)
  maximum = max(phenotype_col, na.rm=T) 
  denominator = maximum - minimum
  print(maximum)
  print(minimum)
  result = sapply(phenotype_col, function(x) ((x - minimum) / denominator))
  print(result)
  return(result)
  
}

code_outliers <- function(vector){
  print(summary(vector))
  maxnonoutlier = 1.5*(summary(vector)[5] - summary(vector)[2]) + summary(vector)[5]
  newvector = if_else(vector > maxnonoutlier, 1, 0)
  newvector[newvector=="NA"] <- 0
  print(newvector)
  return(newvector)
}

phenotypes$hemolysis = sapply( phenotypes$hemolysis, function(x) as.numeric(as.character(x)))


# Max-min normalization
#########################
phenotypes$HemolysisNorm = MinMaxNormalization(phenotypes$hemolysis)
plotheme = ggplot(phenotypes, aes(x=(HemolysisNorm))) + geom_histogram(binwidth=.05,color="black", fill="firebrick4") + xlab("Hemolysis Phenotype ")

phenotypes$BioFilmNorm = MinMaxNormalization(phenotypes$biofilm)
plotfilm = ggplot(phenotypes, aes(x=(BioFilmNorm))) + geom_histogram(binwidth=.05,color="black", fill="darkseagreen") + xlab("Biofilm Phenotype")

phenotypes$XanthinNorm = MinMaxNormalization(phenotypes$xanthin)
plotxan = ggplot(phenotypes, aes(x=(XanthinNorm))) + geom_histogram(binwidth=.05,color="black", fill="gold3") + xlab("Staphyloxanthin Phenotype")
gridExtra::grid.arrange(plotheme, plotfilm, plotxan ,ncol=3)

phenotypes$KinaseNorm = MinMaxNormalization(phenotypes$kinase)
ggplot(phenotypes, aes(x=(KinaseNorm))) + geom_histogram(binwidth=.01,color="black", fill="darkmagenta") + xlab("Staphylokinase Phenotype (Max/Min Normalized)")

phenotypes$SiderophoreNorm = MinMaxNormalization(phenotypes$siderophore)
ggplot(phenotypes, aes(x=(SiderophoreNorm))) + geom_histogram(binwidth=.05,color="black", fill="darkmagenta") + xlab("Siderophore Phenotype (Max/Min Normalized)")


phenotypes$DORNID = paste0("DORN", phenotypes$DORN)

phenotypes$hemolysis_outlier = code_outliers(phenotypes$hemolysis_dev)
phenotypes$biofilm_outlier = code_outliers(phenotypes$biofilm_dev)
phenotypes$xanthin_outlier = code_outliers(phenotypes$xanthin_dev)
phenotypes$kinase_outlier  = code_outliers(phenotypes$kinase_dev)
phenotypes$sum_outliers =  phenotypes$hemolysis_outlier +  phenotypes$biofilm_outlier + phenotypes$xanthin_outlier +  phenotypes$kinase_outlier
View(phenotypes)


# Coefficient of variance 
#########################
phenotypes$hemolysis_coeff_var =  (phenotypes$hemolysis_dev)/mean(phenotypes$hemolysis_dev, na.rm=T)
phenotypes$biofilm_coeff_var =  (phenotypes$biofilm_dev) / mean(phenotypes$biofilm_dev, na.rm=T)
phenotypes$xanthin_coeff_var=  (phenotypes$xanthin_dev) / mean(phenotypes$xanthin_dev, na.rm=T)
phenotypes$kinase_coeff_var =  (phenotypes$kinase_dev) / mean(phenotypes$kinase_dev, na.rm=T)
phenotypes$siderophore_coeff_var =   (phenotypes$siderophore_dev) / mean(phenotypes$siderophore_dev, na.rm=T)

normalized_phenotypes_devs_melted = phenotypes[c("DORN" ,"hemolysis_coeff_var", "biofilm_coeff_var", "xanthin_coeff_var", "kinase_coeff_var", "siderophore_coeff_var")] %>% reshape2::melt(id.vars=c("DORN"))
ggplot(normalized_phenotypes_devs_melted, aes(x=variable, y=value))+ geom_point() 

hemolysisThirdQuartile = (summary(phenotypes$hemolysis_coeff_var))[5]
siderophoreThirdQuartile = (summary(phenotypes$siderophore_coeff_var))[5]
biofilmThirdQuartile = (summary(phenotypes$biofilm_coeff_var))[5]
xanthinThirdQuartile = (summary(phenotypes$xanthin_coeff_var))[5]
kinaseThirdQuartile = (summary(phenotypes$kinase_coeff_var))[5]
phenotypes$sum_outliers

varplot = ggplot(normalized_phenotypes_devs_melted, aes(x=variable, y=value))+ geom_point() +
  geom_segment(aes(x=.55,xend=1.45, y=hemolysisThirdQuartile, yend=hemolysisThirdQuartile), color="red", linetype="dashed")+
  geom_segment(aes(x=1.55,xend=2.45, y=biofilmThirdQuartile, yend=biofilmThirdQuartile), color="red", linetype="dashed")+
  geom_segment(aes(x=2.55,xend=3.45, y=xanthinThirdQuartile, yend=xanthinThirdQuartile), color="red", linetype="dashed") + 
  geom_segment(aes(x=3.55,xend=4.45, y=kinaseThirdQuartile, yend=kinaseThirdQuartile), color="red", linetype="dashed") + 
  geom_segment(aes(x=4.55,xend=5.45, y=siderophoreThirdQuartile, yend=siderophoreThirdQuartile), color="red", linetype="dashed") +
  xlab("Phenotype") + ylab("Coefficient of Variation (SD/Mean)") + ggtitle("Variation in Phenotypes Labeled at 75th Percentile")

  # geom_hline(yintercept=hemolysisThirdQuartile, linetype="solid", color = "red") +
  # geom_hline(yintercept=siderophoreThirdQuartile, linetype="solid", color = "purple") +
  # geom_hline(yintercept=biofilmThirdQuartile, linetype="solid", color = "green") +
  # geom_hline(yintercept=xanthinThirdQuartile, linetype="dashed", color = "yellow") +
  # geom_hline(yintercept=siderophoreThirdQuartile, linetype="dashed", color = "blue") + 
  # ggtitle("Distributions of SD with 3rd Quartile Cutoff Labeled") 

ggsave(varplot, file="Normalized SDs of Phenotypes.png")
ggsave(varplot + geom_text(label=normalized_phenotypes_devs_melted$DORN, hjust=.01), height=40, width=49, file="variability_phenotype_plot.pdf")


ggplot(normalized_phenotypes_devs_melted, aes(x=variable, y=value))+ geom_point() 


View(phenotypes)



phenotypes_devs_melted = phenotypes[c("DORN" ,"hemolysis_dev", "biofilm_dev", "xanthin_dev", "siderophore_dev")] %>% reshape2::melt(id.vars=c("DORN"))

ggplot(phenotypes_devs_melted, aes(x=variable, y=value))+ geom_point() 




hemolysisThirdQuartile = (summary(phenotypes$hemolysis_dev))[5]
siderophoreThirdQuartile = (summary(phenotypes$siderophore_dev))[5]
biofilmThirdQuartile = (summary(phenotypes$biofilm_dev))[5]
xanthinThirdQuartile = (summary(phenotypes$xanthin_dev))[5]
varplot = ggplot(phenotypes_devs_melted, aes(x=variable, y=value))+ geom_point() +
  geom_hline(yintercept=hemolysisThirdQuartile, linetype="solid", color = "red") +
  #geom_hline(yintercept=siderophoreThirdQuartile, linetype="dashed", color = "purple") +
  geom_hline(yintercept=biofilmThirdQuartile, linetype="solid", color = "green") +
  geom_hline(yintercept=hemolysisThirdQuartile, linetype="dashed", color = "yellow") + ggtitle("Distributions of SD with 3rd Quartile Cutoff Labeled") +
  geom_text(label=phenotypes_devs_melted$DORN, hjust=.01)
varplot

ggsave(varplot, height=40, width=49, file="variability_phenotype_plot.pdf")
ggplot(phenotypes, aes(x=biofilm_dev, y=hemolysis_dev))+ geom_point() + geom_text(label=phenotypes$DORNID)
phenotypes$heme_and_biofilm = (phenotypes$hemolysis_outlier==1 & phenotypes$biofilm_outlier==1)
  
ggplot(phenotypes, aes(x=biofilm_dev, y=hemolysis_dev)) + geom_point() + geom_text(label=phenotypes$DORNID)
phenotypes$three = (phenotypes$hemolysis_outlier==1 & phenotypes$biofilm_outlier==1 & phenotypes$kinase==1 )


p1 = ggplot(phenotypes, aes(x=BioFilmNorm, y=XanthinNorm)) + geom_point()  + xlab ("Biofilm Phenotype") + ylab("Staphyloxanthin Phenotype") #+ #geom_text(label=phenotypes$DORNID)
p2 = ggplot(phenotypes, aes(x=HemolysisNorm, y=XanthinNorm)) + geom_point() + xlab ("Hemolysis Phenotype") + ylab("Staphyloxanthin Phenotype")#+ #geom_text(label=phenotypes$DORNID)
p3 =ggplot(phenotypes, aes(x=BioFilmNorm, y=HemolysisNorm)) + geom_point() + xlab ("Biofilm Phenotype") + ylab("Hemolysis Phenotype")#+ #geom_text(label=phenotypes$DORNID)
gridExtra::grid.arrange(p1, p2, p3, ncol=3)



pairs(phenotypes_devs)


ggplot(phenotypes, aes(x=siderophore_dev, y=hemolysis_dev)) + geom_point() + geom_text(label=phenotypes$DORNID)
phenotypes


