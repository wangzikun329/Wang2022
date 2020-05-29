
analysis <- "edgeR_TC_"

# set up
setwd("~/Documents/Laboratory/RNA_seq")
library(edgeR)
# install edgeR from bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")

# read in count table
read.counts <- read.table(file = "Filtered_raw_counts.txt", header = TRUE)
# place WT2U as the first group
read.counts <- read.counts[ , c(c(37:48),c(1:36))]
head(read.counts)

# make design table
sample.info <- data.frame(condition = c(rep("aWT2U", 12), rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12)), 
                          ZT = rep( c(paste("ZT", c("00", "04", "08", "12", "16", "20", "00", "04", "08", "12", "16", "20"),sep="")), 4), 
                          row.names = names(read.counts) )
sample.info
design <- model.matrix(~condition + ZT, data=sample.info)
colnames(design)

# specify sample types
sample.type <- factor(c(rep("aWT2U", 12), rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12)))
sample.type

# convert the count matrix into an edgeR object using DGEList
edgeR.DGElist <- DGEList(counts = read.counts, group = sample.type)

# check the result
head(edgeR.DGElist$counts)
edgeR.DGElist$samples

# calculate normalization factors for the library sizes
edgeR.DGElist <- calcNormFactors( edgeR.DGElist, method = "TMM" )
edgeR.DGElist$samples

# relevel the edgeR object
edgeR.DGElist$samples$group <- relevel(edgeR.DGElist$samples$group, ref="aWT2U")

# generalize linear model (GLM) method

# estimate common dispersion, trended dispersions and tagwise dispersions
edgeR.DGElist <- estimateDisp(edgeR.DGElist, design)

# generalized linear model (GLM) likelihood ratio test
fit <- glmFit(edgeR.DGElist, design)
colnames(fit)

# extract the results, write to file
i <- 2
for (genotype in c("INC1", "WAKED2", "WT2U_old")){
  DGE.results <- glmLRT(fit, coef=i)
  i <- i+1
  out.df <- topTags(DGE.results, n = Inf, sort.by = "PValue", adjust.method = "BH", p.value = 0.05)
  write.table(out.df, file = paste(analysis,genotype,".txt",sep=""), sep = "\t", quote = FALSE, row.names = TRUE)
}
