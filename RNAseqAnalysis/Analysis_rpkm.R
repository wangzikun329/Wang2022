
# setup
setwd("~/Documents/Laboratory/RNA_seq")
library(ggplot2)
library(edgeR)

# read in raw data and info data
read.counts <- read.table(file = "Raw_data_combined.txt", header = TRUE)
read.info <- read.table(file = "Raw_data_combined_info.txt", header = TRUE)

filtered.counts <- read.table(file = "Filtered_sfnorm_counts.txt", header = TRUE)
filtered.info <- read.info[row.names(filtered.counts), ]
write.table(filtered.info, file = "Filtered_counts_info.txt", sep = "\t", quote = FALSE, row.names = TRUE)

counts.rpkm <- rpkm(read.counts, gene.length = read.info$Length, log=FALSE)
filtered.rpkm <- rpkm(filtered.counts, gene.length = filtered.info$Length, log = FALSE)

write.table(counts.rpkm, file = "Raw_data_combined_rpkm.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(filtered.rpkm, file = "Filtered_sfnorm_rpkm.txt", sep = "\t", quote = FALSE, row.names = TRUE)
