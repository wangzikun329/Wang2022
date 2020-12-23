
analysis <- "DESeq2_"

# setup
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ggplot2)
library(DESeq2)
library(edgeR)

# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Dm.eg.db")
# biocLite("DESeq2")
# biocLite("vsn")
# biocLite("NMF")
# biocLite("GO.db")

# read in raw data
read.counts <- read.table(file = "Filtered_raw_counts.txt", header = TRUE)
# place WT2U as the first group
read.counts <- read.counts[ , c(c(37:48),c(1:36))]
names(read.counts)

# define sample metadata
sample.info <- data.frame(condition = c(rep("aWT2U", 12), rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12)), row.names = names(read.counts) )
sample.info

# generate DESeqDataSet
dseq.ds <- DESeqDataSetFromMatrix(count = read.counts, colData = sample.info, design = ~ condition)
head(counts(dseq.ds))
colSums(counts(dseq.ds))

# relevel the dseq.ds dataframe to compare to WT2U
dseq.ds$condition <- relevel(dseq.ds$condition, "aWT2U")
dseq.ds$condition

# run DESeq
dseq.ds <- DESeq(dseq.ds)

# extract and save DESeq result
for (genotype in c("INC1", "WAKED2", "WT2U_old")){
  DGE.results <- results(dseq.ds, independentFiltering = TRUE, pAdjustMethod = "BH", contrast = c("condition", genotype, "aWT2U"))
  # summary(DGE.results)
  # table(DGE.results$padj < 0.05)
  out.df <- subset( as.data.frame(DGE.results), padj < 0.05)
  write.table(out.df, file = paste(analysis,genotype,".txt",sep=""), sep = "\t", quote = FALSE, row.names = TRUE)
}
