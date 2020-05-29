
analysis <- "DESeq2_TC_"

# setup
setwd("~/Documents/Laboratory/RNA_seq")
library(ggplot2)
library(DESeq2)

# read in raw data
read.counts <- read.table(file = "Filtered_raw_counts.txt", header = TRUE)
# place WT2U as the first group
read.counts <- read.counts[ , c(c(37:48),c(1:36))]
names(read.counts)

# make design table
sample.info <- data.frame(condition = c(rep("aWT2U", 12), rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12)), 
                          ZT = rep( c(paste("ZT", c("00", "04", "08", "12", "16", "20", "00", "04", "08", "12", "16", "20"),sep="")), 4), 
                          row.names = names(read.counts) )
sample.info

# generate DESeqDataSet
dseq.ds <- DESeqDataSetFromMatrix(count = read.counts, colData = sample.info, design = ~ condition + ZT)

head(counts(dseq.ds))
colSums(counts(dseq.ds))

# calculate and correct for differences in library sizes
# dseq.ds <- estimateSizeFactors(dseq.ds)
# sizeFactors(dseq.ds)

# relevel the dseq.ds dataframe to compare to WT2U
dseq.ds$condition <- relevel(dseq.ds$condition, "aWT2U")
dseq.ds$condition

# run DESeq
dseq.ds <- DESeq(dseq.ds, test="LRT", reduced = ~ ZT)
resultsNames(dseq.ds)

# extract and save DESeq result
for (genotype in c("INC1", "WAKED2", "WT2U_old")){
  DGE.results <- results(dseq.ds, independentFiltering = TRUE, pAdjustMethod = "BH", name = paste("condition",genotype,"vs_aWT2U",sep="_"), test = "Wald")
  # summary(DGE.results)
  # table(DGE.results$padj < 0.05)
  out.df <- subset( as.data.frame(DGE.results), padj < 0.05)
  write.table(out.df, file = paste(analysis,genotype,".txt",sep=""), sep = "\t", quote = FALSE, row.names = TRUE)
}

# single gene plot with ZT as x
gene.data <- plotCounts(dseq.ds, gene = "FBgn0023076", normalized = TRUE, intgroup=c("condition","ZT"), returnData=TRUE)
ggplot(gene.data, aes(x=ZT, y=count, color=condition, group=condition)) + geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()



