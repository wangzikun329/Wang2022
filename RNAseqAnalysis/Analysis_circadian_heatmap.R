
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(NMF)

# read in circadian count data and phase data
circadian.data <- read.table(file = "circadian_union_WT2U_result.txt", header = TRUE)
jtk.out <- read.table(file = "JTK_CYCLE_WT2U_out.txt", header = TRUE)
row.names(jtk.out) <- jtk.out$x
read.rpkm <- read.table(file = "Raw_data_combined_rpkm.txt", header = TRUE)

# extract subset of information from jtk.out
phase <- jtk.out[row.names(circadian.data), ]

# order phase and data according to phase
phase <- phase[order(phase$LAG),]
circadian.data <- read.rpkm[row.names(phase),]

# draw heatmap for each genotype
genotype <- c(rep("INC1", 12), rep("WAKED2", 12), rep("WT2U_old", 12), rep("WT2U", 12))
for (i in c(1,13,25,37)){
  plot.data <- circadian.data[, i:(i+11)]
  aheatmap(plot.data, Rowv = NA, Colv = NA, scale = "row")
  aheatmap(plot.data, Rowv = NA, Colv = NA, scale = "row", filename = paste("heatmap_all_circadian_", genotype[i], ".pdf", sep=""))
}

# save wt2u plot
wt2u.data <- circadian.data[, 37:48]
aheatmap(wt2u.data, Rowv = NA, Colv = NA, scale = "row", filename = "heatmap_all_circadian.pdf")

# save ordered circadian data with JTK result
out.df <- cbind(phase[ , 2:6], circadian.data)
# write.table(out.df, file = "circadian_union_JTK_WT2U.txt", sep = "\t", quote = FALSE, row.names = TRUE)
