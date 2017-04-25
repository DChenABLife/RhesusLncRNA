#########################################################################################
#                                                                                       #
#   Copyright (c)   A_B_Life 2014                                                       #
#   Writer:         Cheng Chao （chaocheng@ablife.cc）                                  #
#   Program Date:   2014.05.19                                                          #
#   Modifier:       Cheng Chao （chaocheng@ablife.cc）                                  #
#   Last Modified:  2014.05.19                                                          #
#                                                                                       #
#   if have some questions,contact the Writer or Modifier                               #
#                                                                                       #
#########################################################################################

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop('Your input arguments is wrong!\n
		args1:\t the mRNA expression profile\n
		args2:\t the miRNA expression profile\n
		args3:\t the correlation result file\n
		args4:\t the pvalue result file\n
		args5:\t the correlation method,pearson:0,spearman:1\n
		args6:\t the out directory\n\n')
}
mrnafile <- args[1]
mirnafile <- args[2]
outfile <- args[3]
outfile2 <- args[4]
cor_method <- as.numeric(args[5])
outdir <- args[6]

setwd(outdir)

library(psych)

# mrnadata <- read.table(file=mrnafile,header=T,sep='\t',row.names="Gene")
# mirnadata <- read.table(file=mirnafile,header=T,sep='\t',row.names="miRNA")

mrnadata <- read.table(file=mrnafile,header=T,sep='\t',row.names=1)
mirnadata <- read.table(file=mirnafile,header=T,sep='\t',row.names=1)

if (cor_method == 0) {
	# result <- corr.test(t(mrnadata),t(mirnadata), use = "complete",method="pearson",adjust="none")
	result <- corr.test(t(mrnadata),t(mirnadata), use = "complete",method="pearson",adjust="none", ci=FALSE)
}
if (cor_method == 1) {
	result <- corr.test(t(mrnadata),t(mirnadata), use = "complete",method="spearman",adjust="holm", ci=FALSE)
}


colnames(result$p)[1] <- paste('mRNA/miRNA(pvalue)',colnames(result$p)[1],sep="\t")
colnames(result$r)[1] <- paste('mRNA/miRNA(correlation)',colnames(result$r)[1],sep="\t")
# colnames(result$p)[1]
# colnames(result$r)[1]

write.table(result$r,file=outfile,row.names=T,col.names=TRUE,append=F,quote=F,sep='\t')
write.table(result$p,file=outfile2,row.names=T,col.names=TRUE,append=F,quote=F,sep='\t')
