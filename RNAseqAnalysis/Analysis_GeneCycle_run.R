
setwd("~/Documents/Laboratory/RNA_seq")

project <- "WT2U"

# install GeneCycle package
# install.packages("fdrtool_1.2.15.tar.gz",repos = NULL, type = "source")
# install.packages("corpcor_1.6.8.tar.gz",repos = NULL, type = "source")
# install.packages("longitudinal_1.1.12.tar.gz",repos = NULL, type = "source")
# install.packages("GeneCycle_1.1.2.tar.gz",repos = NULL, type = "source")

library("GeneCycle")

# read in count data
data <- read.table(file = paste("JTK_CYCLE_",project,".txt",sep=""), header = TRUE)

# save gene names
anno <- row.names(data)

# transpose data
data.input <- t(data)

# Fishers exact G test
pval.data <- fisher.g.test(data.input)

# extract the results
fdr.out <- fdrtool(pval.data, statistic = "pvalue")

# put together total results
result.total <- cbind(as.data.frame(fdr.out), data)
out.df <- subset( as.data.frame(result.total), pval < 0.05)

# save result
write.table(result.total,file=paste("GeneCycle_",project,"_out.txt",sep=""),row.names=T,col.names=T,quote=F,sep="\t")
write.table(out.df,file=paste("GeneCycle_",project,"_result.txt",sep=""),row.names=T,col.names=T,quote=F,sep="\t")

