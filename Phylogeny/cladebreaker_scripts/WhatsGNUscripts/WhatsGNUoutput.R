# Amy Campbell
# 3/2023
# reading output from WhatsGNU plotter function

args <- commandArgs(trailingOnly = TRUE)

PathToWhatsGNU=args[1]

WhatsGNUvolcano=read.csv2(PathToWhatsGNU, sep="\t")

outputpath = dirname(PathToWhatsGNU)

outputfile = paste0(outputpath,"/", basename(outputpath), "_ortholog_list.txt")


WhatsGNUvolcano$X.log10.pvalue = sapply(WhatsGNUvolcano$X.log10.pvalue, function(x) as.numeric(as.character(x)))
WhatsGNUvolcano$DELTA_AVG_GNU = sapply(WhatsGNUvolcano$DELTA_AVG_GNU, function(x) as.numeric(as.character(x)))
TopOrthologs= WhatsGNUvolcano %>% filter((X.log10.pvalue >2)) %>% filter(abs(DELTA_AVG_GNU)>5000) %>% filter(case==max(case) & control==max(control))
write.table(TopOrthologs$ortho_group, quote = F, row.names=F, col.names=F, file=outputfile)
