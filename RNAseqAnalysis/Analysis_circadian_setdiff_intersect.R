
setwd("~/Documents/Laboratory/RNA_seq")
library(NMF)

# read in circadian count data and phase data
inc1.data <- read.table(file = "circadian_setdiff_INC1_out.txt", header = TRUE)
waked2.data <- read.table(file = "circadian_setdiff_WAKED2_out.txt", header = TRUE)
wt2uold.data <- read.table(file = "circadian_setdiff_WT2U_old_out.txt", header = TRUE)

inc1.loss <- inc1.data[1:161,]
waked2.loss <- waked2.data[1:153,]
wt2uold.loss <- wt2uold.data[1:111,]

inc1.genename <- row.names(inc1.loss)
waked2.genename <- row.names(waked2.loss)
wt2uold.genename <- row.names(wt2uold.loss)

intersect.genename <- intersect(intersect(inc1.genename, waked2.genename), wt2uold.genename)

# plot all the genes when dseq.ds is already available
# for (genename in intersect.genename){
#   rpkmPlotbyID(genename)
# }
# 

read.rpkm <- read.table(file = "Raw_data_combined_rpkm.txt", header = TRUE)
out.data <- read.rpkm[intersect.genename,]
write.table(out.data, file = "circadian_setdiff_RhythmicityLoss.txt", quote = FALSE, row.names = TRUE)
