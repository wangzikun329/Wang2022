
# setup
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ggplot2)
library(reshape2)
library(DESeq2)
library(edgeR)

# read in raw data
read.counts <- read.table(file = "Raw_data_combined.txt", header = TRUE)
names(read.counts)

# define sample metadata
sample.info <- data.frame(condition = c(rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12), rep("WT2U", 12)), 
                          ZT = rep( c(paste("ZT", c("00", "04", "08", "12", "16", "20", "00", "04", "08", "12", "16", "20"),sep="")), 4), 
                          row.names = names(read.counts) )
sample.info

# generate DESeqDataSet
dseq.ds <- DESeqDataSetFromMatrix(count = read.counts, colData = sample.info, design = ~ condition + ZT)
head(counts(dseq.ds))
colSums(counts(dseq.ds))

# remove genes without any reads in any condition
# dseq.ds <- dseq.ds[rowSums(counts(dseq.ds)) > 0,]

# find the threshold for low expressed genes
for(j in seq(0.1, 1, length = 10)){
  low.expressed <- data.frame(row.sum = 1:48, keep.count = rep(0,48))
  for(i in low.expressed[,"row.sum"]){
    keep <- rowSums( cpm(counts(dseq.ds)) >= j) >= i
    low.expressed[i, "keep.count"] <- count(keep)
  }
  plot(low.expressed)
  title(main = paste("row sums of cpm counts >= ", j))
}

# remove genes with low expression level
keep <- rowSums( cpm(counts(dseq.ds)) >= 1) >= 20
dseq.ds <- dseq.ds[keep, ]
  

# calculate and correct for differences in library sizes
dseq.ds <- estimateSizeFactors(dseq.ds)
sizeFactors(dseq.ds)

# extract raw counts and size factor normalized counts
counts.raw <- data.frame(counts(dseq.ds))
counts.sfnorm <- data.frame(counts(dseq.ds, normalized = TRUE))

# save fileted data for next step analysis
# write.table(counts.raw, file = "Filtered_raw_counts.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# write.table(counts.sfnorm, file = "Filtered_sfnorm_counts.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# check the similarities between samples
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

# function to heatmap correlation between samples
heatmap.correlation <- function(data, digits = 4, limit = c(0.8,1), text_label = FALSE)
{
  # calculate the Pearson Correlation between samples
  data.cormat <- round(cor(data), digits = digits)
  melted.cormat <- melt(data.cormat)
  
  # generate heatmap
  ggheatmap <- ggplot(data = melted.cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = mean(limit), limit = limit, space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    scale_y_discrete(limits = rev(levels(melted.cormat[,2]))) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0), axis.title=element_blank())+
    coord_fixed()
  
  # label heatmap palette with correlation 
  if(text_label == TRUE){
    ggheatmap <- ggheatmap +
      geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
  }
  
  print(ggheatmap)
}


# transform size factor normalized counts to log2 scale
log.norm.counts <- log2(counts.sfnorm + 0.001)

# plot the similarities within groups
# pairs(log.norm.counts[,c(paste("INC1",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.smooth, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(log.norm.counts[,c(paste("WAKED2",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.smooth, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(log.norm.counts[,c(paste("WT2U_old",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.smooth, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(log.norm.counts[,c(paste("WT2U",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.smooth, upper.panel = panel.cor, xaxt = "n", yaxt = "n")

pairs(log.norm.counts[,c(paste("INC1",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
pairs(log.norm.counts[,c(paste("WAKED2",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
pairs(log.norm.counts[,c(paste("WT2U_old",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
pairs(log.norm.counts[,c(paste("WT2U",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")

# WAKED2_ZT12 and WAKED2_ZT36 are the two replicates with the lowest similarities
pairs(log.norm.counts[,c("WAKED2_ZT12", "WAKED2_ZT36")], lower.panel = panel.smooth, upper.panel = panel.cor)

# the rlog transformation will decrease the variation for lowly expressed genes.
rlog.dseq.sumExp <- rlog(dseq.ds, blind = TRUE)

pairs(assay(rlog.dseq.sumExp[,c(paste("INC1",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
pairs(assay(rlog.dseq.sumExp[,c(paste("WAKED2",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
pairs(assay(rlog.dseq.sumExp[,c(paste("WT2U_old",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
pairs(assay(rlog.dseq.sumExp[,c(paste("WT2U",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")

# vst is an alternative and faster estimation of the dispersion trend
# vst.dseq.sumExp <- vst(dseq.ds, blind = TRUE, nsub = 1000, fitType = "parametric")
# 
# pairs(assay(vst.dseq.sumExp[,c(paste("INC1",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(assay(vst.dseq.sumExp[,c(paste("WAKED2",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(assay(vst.dseq.sumExp[,c(paste("WT2U_old",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(assay(vst.dseq.sumExp[,c(paste("WT2U",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))]), lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")

# plot transformed vs nontransformed using the most different replicates
plot(assay(rlog.dseq.sumExp[,16]), assay(rlog.dseq.sumExp[,22]))
# plot(assay(vst.dseq.sumExp[,16]), assay(vst.dseq.sumExp[,22]))
plot(log.norm.counts[,16], log.norm.counts[,22])

# heatmap correlation between samples
heatmap.correlation(log.norm.counts, limit = c(0.8,1))
heatmap.correlation(assay(rlog.dseq.sumExp), limit = c(0.95,1))
# heatmap.correlation(assay(vst.dseq.sumExp), limit = c(0.9,1))

plotPCA_with_label <- function(dseq.sum, ntop = 500, returnData = FALSE){
  pca <- plotPCA(dseq.sum, ntop = ntop, returnData = TRUE)
    ZT <- data.frame(ZT = rep( c(paste("ZT", c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="")), 4),
                     row.names = row.names(pca))
  pca <- cbind(pca, ZT)
  pca.plot <- ggplot(pca, aes(x=PC1, y=PC2, color=condition, group=condition, label = ZT)) + 
    geom_point() + geom_text(hjust = 0, nudge_x = 0.5)
  print(pca.plot)
  if(returnData == TRUE){
    return(pca)
  }
}

# plot PCA
pca.results <- plotPCA_with_label(rlog.dseq.sumExp, 2000, returnData = TRUE)
# plotPCA(vst.dseq.sumExp, ntop = 500)
# save PCA results
write.table(pca.results, file = "PCA_rlog_transfromed_results.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# inspect the within-group pca
plotPCA_single <- function(dseq.sum, ntop = 1000, returnData = FALSE){
  pca <- plotPCA(dseq.sum, ntop = ntop, returnData = TRUE)
  ZT <- data.frame(ZT = c(paste("ZT", c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="")),
                   row.names = row.names(pca))
  pca <- cbind(pca, ZT)
  pca.plot <- ggplot(pca, aes(x=PC1, y=PC2, color=name, group=name, label = ZT)) + 
    geom_point() + geom_text(hjust = 0, nudge_x = 0.1)
  print(pca.plot)
  if(returnData == TRUE){
    return(pca)
  }
}
plotPCA_single(rlog.dseq.sumExp[,1:12])
plotPCA_single(rlog.dseq.sumExp[,13:24])
plotPCA_single(rlog.dseq.sumExp[,25:36])
plotPCA_single(rlog.dseq.sumExp[,37:48])

# the code below takes into consideration that some genes in wakeD2 showed abnormal expression
# these genes showed higher expression in a handful of timepoints irrespective to 24h cycle
# these genes will be removed
# calculate differences between replicates
# diff.std.error <- counts.sfnorm[,13:24] / (rowSums(counts.sfnorm[,13:24])/12)
# diff.std.error <- log2(diff.std.error + 0.001)
# gene.list <- log.norm.counts[abs(rowSums(diff.std.error)) >= 10,]
# genes.large.variance <- intersect(row.names(gene.list), row.names(wake.zt36))


# max.diffrence <- log.norm.counts[,16]-log.norm.counts[,22]
# hist(max.diffrence, breaks = 200)
# sum(abs(max.diffrence) > 10)
# 
# max.difference.thres <- 1:20
# for(i in 1:20){
#   max.difference.thres[i] <- sum(abs(max.diffrence) > i)
# }
# plot(max.difference.thres)
# 
# # save genes to be discarded
# genes.large.variance <- counts.raw[abs(max.diffrence) > 5, ]
# write.table(genes.large.variance, file = "Genes_filtered_by_variance.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# 
# for(gene in intersect(row.names(genes.large.variance), row.names(dge.tc))) rpkmPlotbyID(gene)
# 
# for(gene in setdiff(row.names(dge.tc), row.names(counts.raw.filter))) rpkmPlotbyID(gene)
# 
# # filter normalized counts with threshold of difference
# log.norm.filter <- log.norm.counts[abs(max.diffrence) <= 5,]
# 
# # apply same filter to deseq data set
# dseq.ds.filter <- dseq.ds[row.names(log.norm.filter), ]
# 
# # calculate and correct for differences in library sizes
# dseq.ds.filter <- estimateSizeFactors(dseq.ds.filter)
# sizeFactors(dseq.ds.filter)
# 
# # extract raw counts and size factor normalized counts
# counts.raw.filter <- data.frame(counts(dseq.ds.filter))
# counts.sfnorm.filter <- data.frame(counts(dseq.ds.filter, normalized = TRUE))
# 
# # transform size factor normalized counts to log2 scale
# log.norm.counts.filter <- log2(counts.sfnorm.filter + 0.001)
# 
# # rlog transform filtered deseq data set
# rlog.dseq.filter <- rlog(dseq.ds.filter, blind = TRUE)
# vst.dseq.filter <- vst(dseq.ds.filter, blind = TRUE, nsub = 1000, fitType = "parametric")
# 
# # plot rlog transformed vs nontransformed
# plot(assay(rlog.dseq.filter[,16]), assay(rlog.dseq.filter[,22]))
# plot(assay(vst.dseq.filter[,16]), assay(vst.dseq.filter[,22]))
# plot(log.norm.counts.filter[,16], log.norm.counts.filter[,22])
# 
# # plot PCA
# plotPCA(rlog.dseq.filter)
# plotPCA(vst.dseq.filter)
# 
# heatmap.correlation(log.norm.filter, limit = c(0.8,1))
# heatmap.correlation(assay(rlog.dseq.filter), limit = c(0.95,1))
# heatmap.correlation(assay(vst.dseq.filter), limit = c(0.9,1))
# 
# plotPCA(rlog.dseq.filter[,1:12])
# plotPCA(rlog.dseq.filter[,13:24])
# plotPCA(rlog.dseq.filter[,25:36])
# plotPCA(rlog.dseq.filter[,37:48])
# 
# plotPCA(vst.dseq.filter[,1:12])
# plotPCA(vst.dseq.filter[,13:24])
# plotPCA(vst.dseq.filter[,25:36])
# plotPCA(vst.dseq.filter[,37:48])
# 
# # plot the similarities within groups after filterting
# pairs(log.norm.filter[,c(paste("INC1",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(log.norm.filter[,c(paste("WAKED2",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(log.norm.filter[,c(paste("WT2U_old",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# pairs(log.norm.filter[,c(paste("WT2U",c("00", "04", "08", "12", "16", "20", "24", "28", "32", "36", "40", "44"),sep="_ZT"))], lower.panel = panel.cor, upper.panel = panel.cor, xaxt = "n", yaxt = "n")
# 
# write.table(counts.raw.filter, file = "Filtered_raw_counts_PCA.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# write.table(counts.sfnorm.filter, file = "Filtered_sfnorm_counts_PCA.txt", sep = "\t", quote = FALSE, row.names = TRUE)

