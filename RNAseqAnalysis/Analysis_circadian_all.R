
project <- "WT2U"

# setup
setwd("~/Documents/Laboratory/RNA_seq")

# read in results
arser.out <- read.table(file = paste("ARSER_", project, "_out.txt", sep=""),header = TRUE)
row.names(arser.out) <- arser.out$probe
jtk.out <- read.table(file = paste("JTK_CYCLE_", project, "_out.txt", sep=""), header = TRUE)
row.names(jtk.out) <- jtk.out$x

head(arser.out)
head(jtk.out)

# subset result to save significant data
out.df <- subset( as.data.frame(arser.out), pvalue < 0.05)

# write.table(out.df,file=paste("ARSER_", project, "_result.txt", sep=""),row.names=F,col.names=T,quote=F,sep="\t")

arser.result <- read.table(file = paste("ARSER_", project, "_result.txt", sep=""),header = TRUE)

jtk.result <- read.table(file = paste("JTK_CYCLE_",  project, "_result.txt", sep=""),header = TRUE)

genecycle.result <- read.table(file = paste("GeneCycle_",  project, "_result.txt", sep=""),header = TRUE)

head(jtk.result)

intersect.result <- intersect(intersect(arser.result$probe, jtk.result$x), row.names(genecycle.result))

arser.jtk <- intersect(arser.result$probe, jtk.result$x)
arser.genecycle <- intersect(arser.result$probe, row.names(genecycle.result))
jtk.genecycle <- intersect(jtk.result$x, row.names(genecycle.result))

temp <- intersect(jtk.genecycle, arser.jtk)

union.result <- union(union(arser.genecycle, arser.jtk), jtk.genecycle)

read.counts <- read.table(file = "Filtered_sfnorm_counts.txt", header = TRUE)

out.union <- read.counts[union.result, ]
out.intersect <- read.counts[intersect.result, ]

# write.table(out.union, file = paste("circadian_union_", project, "_result.txt", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)
# write.table(out.intersect, file = paste("circadian_intersect_", project, "_result.txt", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)

phase <- jtk.out[union.result, ]
union.jtk <- cbind(phase[ , 2:6], out.union)
write.table(union.jtk, file = paste("circadian_union_JTK_", project, "_result.txt", sep = ""), sep = "\t", quote = FALSE, row.names = TRUE)
