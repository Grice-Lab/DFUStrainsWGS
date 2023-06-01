
data531 = read.csv("Documents/DataInputGithub/data/Validations/AmySTX5-31.csv")
ggplot(data531, aes(x=factor(Strain), y=log(PctJE2))) + geom_boxplot() + ggpubr::stat_compare_means(method="t.test", ref.group="632")

