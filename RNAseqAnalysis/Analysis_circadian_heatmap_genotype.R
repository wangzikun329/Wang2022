
setwd("~/Documents/Laboratory/RNA_seq")
library(NMF)

# read in circadian count data and phase data
wt2u.union <- read.table(file = "circadian_union_JTK_WT2U_result.txt", header = TRUE)
inc1.union <- read.table(file = "circadian_union_JTK_INC1_result.txt", header = TRUE)
waked2.union <- read.table(file = "circadian_union_JTK_WAKED2_result.txt", header = TRUE)
wt2uold.union <- read.table(file = "circadian_union_JTK_WT2U_old_result.txt", header = TRUE)

wt2u.data <- wt2u.union[order(wt2u.union$LAG), 42:53]
inc1.data <- inc1.union[order(inc1.union$LAG), 6:17]
waked2.data <- waked2.union[order(waked2.union$LAG), 18:29]
wt2uold.data <- wt2uold.union[order(wt2uold.union$LAG), 30:41]

aheatmap(wt2u.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 0.38, cellwidth = 16, 
         annColors = "RdBu", filename = "heatmap_circadian_WT2U_1115.pdf")
aheatmap(inc1.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 0.38, cellwidth = 16, 
         annColors = "RdBu", filename = "heatmap_circadian_INC1_663.pdf")
aheatmap(waked2.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 0.38, cellwidth = 16, 
         annColors = "RdBu", filename = "heatmap_circadian_WAKED2_809.pdf")
aheatmap(wt2uold.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 0.38, cellwidth = 16, 
         annColors = "RdBu", filename = "heatmap_circadian_WT2U_old_1583.pdf")
