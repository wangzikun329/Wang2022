
project <- "WT2U_old"

setwd("~/Documents/Laboratory/RNA_seq")
library(NMF)

wt2u.union <- read.table(file = "JTK_CYCLE_WT2U_out.txt", header = TRUE)
proj.union <- read.table(file = paste("JTK_CYCLE_",project,"_out.txt",sep=""), header = TRUE)

row.names(wt2u.union) <- wt2u.union$x
row.names(proj.union) <- proj.union$x

# wt2u.zscore <- data.frame(t(scale(t(wt2u.union[7:18]))))
# proj.zscore <- data.frame(t(scale(t(proj.union[7:18]))))
# total.union <- cbind(wt2u.zscore, proj.zscore)

wt2u.high <- wt2u.union[wt2u.union$ADJ.P<=0.01,]
wt2u.ar <- wt2u.union[wt2u.union$ADJ.P>0.5,]
proj.high <- proj.union[proj.union$ADJ.P<=0.01,]
proj.ar <- proj.union[proj.union$ADJ.P>0.5,]

wt2u.only <- intersect(row.names(wt2u.high), row.names(proj.ar))
proj.only <- intersect(row.names(wt2u.ar), row.names(proj.high))

wt2u.only.data <- wt2u.union[wt2u.only,]
wt2u.only.data <- wt2u.only.data[order(wt2u.only.data$LAG),]
proj.only.data <- proj.union[proj.only,]
proj.only.data <- proj.only.data[order(proj.only.data$LAG),]

aheatmap(wt2u.only.data[7:18], Rowv = NA, Colv = NA, scale = "row")
aheatmap(proj.only.data[7:18], Rowv = NA, Colv = NA, scale = "row")

read.rpkm <- read.table(file = "Raw_data_combined_rpkm.txt", header = TRUE)

if(project == "INC1") column.range = c(c(37:48), c(1:12))
if(project == "WAKED2") column.range = c(c(37:48), c(13:24))
if(project == "WT2U_old") column.range = c(c(37:48), c(25:36))

top.data <- read.rpkm[row.names(wt2u.only.data),column.range]
bottom.data <- read.rpkm[row.names(proj.only.data),column.range]

heatmap.data <- rbind(top.data,bottom.data)
row.anno <- c(rep(0, length(wt2u.only)), rep(1, length(proj.only)))

aheatmap(top.data, Rowv = NA, Colv = NA, scale = "row")
aheatmap(bottom.data, Rowv = NA, Colv = NA, scale = "row")
aheatmap(heatmap.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 2, cellwidth = 10, 
         annRow = row.anno, annColors = "RdBu")

# aheatmap(heatmap.data, Rowv = NA, Colv = NA, scale = "row", cellheight = 2, cellwidth = 10, 
#          annRow = row.anno, annColors = "RdBu", filename = paste("circadian_setdiff_",project,".pdf",sep = ""))

write.table(heatmap.data, file = paste("circadian_setdiff",project,"out.txt",sep="_"), sep = "\t", quote = FALSE, row.names = TRUE)
         