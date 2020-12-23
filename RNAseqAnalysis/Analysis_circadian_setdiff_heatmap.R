
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(NMF)
library(org.Dm.eg.db)

# read in circadian count data and phase data
circadian.data <- read.table(file = "circadian_setdiff_RhythmicityLoss.txt", header = TRUE)
circadian.data <- circadian.data[,c(37:48,1:36)]
jtk.out <- read.table(file = "JTK_CYCLE_WT2U_out.txt", header = TRUE)
row.names(jtk.out) <- jtk.out$x
read.rpkm <- read.table(file = "Raw_data_combined_rpkm.txt", header = TRUE)

# extract subset of information from jtk.out
phase <- jtk.out[row.names(circadian.data), ]

# order phase and data according to phase
phase <- phase[order(phase$LAG),]
circadian.data <- circadian.data[row.names(phase),]

anno.data <- read.table(file = "Raw_data_combined_id_curated.txt", header = TRUE, row.names = 1)
anno.data <- anno.data[row.names(circadian.data),]
row.names(circadian.data) <- anno.data$SYMBOL
# wt.mean <- rowMeans(circadian.data[,c(1:12)])
# circadian.data.norm <- circadian.data/wt.mean

aheatmap(circadian.data, Rowv = NA, Colv = NA, scale = "row")
aheatmap(circadian.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 8, cellwidth = 10, 
         annColors = "RdBu", filename = "circadian_setdiff_RhythmmicityLoss.pdf")
